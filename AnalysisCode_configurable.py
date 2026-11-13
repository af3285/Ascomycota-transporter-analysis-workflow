
from Bio import SeqIO, AlignIO, Entrez, pairwise2
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Blast import NCBIXML
import os
import subprocess
from collections import Counter
import argparse
from pathlib import Path
#This pipeline was created by Aaron Fryzel on 6/23/25
#Last updated on 6/25/25
########################################################################################
def extract_and_fetch_sequences(blast_xml, output_fasta):
    # Extract unique protein accessions from all BLAST query records.
    Entrez.email = "af3285@drexel.edu"
    ids = []

    with open(blast_xml) as f:
        for blast_record in NCBIXML.parse(f):
            for hit in blast_record.alignments:
                # Prefer the accession supplied by Biopython when available.
                accession = getattr(hit, "accession", None)

                if not accession:
                    parts = hit.hit_id.split("|")

                    if len(parts) >= 4 and parts[0].lower() == "gi":
                        accession = parts[3]
                    elif len(parts) >= 2 and parts[1]:
                        accession = parts[1]
                    else:
                        accession = hit.hit_id

                if accession and accession not in ids:
                    ids.append(accession)

    # Download sequences from NCBI in batches of 100.
    with open(output_fasta, "w") as out_f:
        for i in range(0, len(ids), 100):
            batch = ids[i:i + 100]
            handle = Entrez.efetch(
                db="protein",
                id=batch,
                rettype="fasta",
                retmode="text"
            )
            try:
                out_f.write(handle.read())
            finally:
                handle.close()

################################################################################ 
#Working with the PDB that was provided           
def extract_sequence_from_pdb(pdb_file, output_fasta):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("query", pdb_file)
    ppb = PPBuilder()
    peptides = ppb.build_peptides(structure)
    if not peptides:
        raise ValueError(f"No polypeptides found in {pdb_file}.")
    #Takes the longest peptide chain
    longest = max(peptides, key=len)
    seq = longest.get_sequence()
    record = SeqRecord(seq, id="query_from_pdb", description="Extracted from PDB")
    SeqIO.write(record, output_fasta, "fasta")


def is_good_alignment(query_seq, target_seq, min_coverage=0.9):
    if not query_seq or not target_seq:
        return False

    # Local protein alignment identifies the region shared by query and target.
    alignments = pairwise2.align.localms(
        query_seq,
        target_seq,
        2,      # Match score
        -1,     # Mismatch score
        -10,    # Gap opening penalty
        -0.5,   # Gap extension penalty
        one_alignment_only=True
    )

    if not alignments:
        return False

    best = alignments[0]

    # Count query residues aligned to an actual target residue.
    covered = sum(
        q != '-' and t != '-'
        for q, t in zip(best.seqA, best.seqB)
    )

    return covered / len(query_seq) >= min_coverage


def filter_by_alignment_coverage(query_fasta, input_fasta, output_fasta, min_coverage=0.9):
    query_record = next(SeqIO.parse(query_fasta, "fasta"))
    query_seq = str(query_record.seq)
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if is_good_alignment(query_seq, str(record.seq), min_coverage):
                SeqIO.write(record, out_f, "fasta")
#################################################################################################################
def filter_sequences(input_fasta, output_fasta, min_length=300):
    # Remove short sequences and sequences containing unknown residues.
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            sequence = str(record.seq).upper()

            if len(sequence) >= min_length and 'X' not in sequence:
                SeqIO.write(record, out_f, "fasta")


def run_cdhit(input_fasta, output_prefix, identity=0.8, coverage=0.9):
    # Cluster sequences at the specified identity and shorter-sequence coverage thresholds.
    subprocess.run([
        "cd-hit",
        "-i", input_fasta,
        "-o", output_prefix,
        "-c", str(identity),
        "-aS", str(coverage),
        "-G", "1",
        "-d", "0",
        "-M", "16000",
        "-T", "4"
    ], check=True)


def parse_clusters(clstr_file, fasta_file, output_file):
    # Use the representative sequence selected by CD-HIT (marked with *).
    representatives = []

    with open(clstr_file) as f:
        for line in f:
            if line.startswith(">Cluster"):
                continue

            if line.rstrip().endswith("*"):
                if ">" not in line:
                    raise ValueError(
                        f"Malformed CD-HIT cluster entry: {line.rstrip()}"
                    )

                description = line.split(">", 1)[1]
                seq_id = description.split("...", 1)[0].split()[0]
                representatives.append(seq_id)

    records = {
        record.id: record
        for record in SeqIO.parse(fasta_file, "fasta")
    }

    missing = [
        seq_id for seq_id in representatives
        if seq_id not in records
    ]

    if missing:
        raise ValueError(
            f"{len(missing)} CD-HIT representative sequence(s) "
            "were not found in the input FASTA file."
        )

    with open(output_file, 'w') as out_f:
        for seq_id in representatives:
            SeqIO.write(records[seq_id], out_f, "fasta")


def run_mafft(input_fasta, output_fasta):
    #Running MAFFT to align the sequences
    with open(output_fasta, 'w') as out_f:
        subprocess.run([
            "mafft", "--maxiterate", "1000", "--localpair", input_fasta # maxiterate == 1000 means that aligments should be high-quality
        ], stdout=out_f, check=True)

def filter_alignment(input_aln, output_aln, max_gap_fraction=0.2):
    # Filter alignment columns that contain too many gaps.
    alignment = AlignIO.read(input_aln, "fasta")

    if len(alignment) == 0:
        raise ValueError("The input alignment contains no sequences.")

    keep_cols = []

    for i in range(alignment.get_alignment_length()):
        col = [record.seq[i] for record in alignment]

        if col.count('-') / len(col) <= max_gap_fraction:
            keep_cols.append(i)

    # Trim each sequence while preserving its record metadata.
    trimmed_records = []

    for record in alignment:
        new_seq = ''.join(str(record.seq[i]) for i in keep_cols)
        record.seq = Seq(new_seq)
        trimmed_records.append(record)

    AlignIO.write(
        MultipleSeqAlignment(trimmed_records),
        output_aln,
        "fasta"
    )


def next_best_consensus(alignment_file, output_file, max_gap_fraction=0.2):
    #Building a consensus sequence from the alignment (skips columns with too many gaps)
    alignment = AlignIO.read(alignment_file, "fasta")
    aln_len = alignment.get_alignment_length()
    consensus = ""
    for i in range(aln_len):
        col = [record.seq[i] for record in alignment]
        if col.count('-') / len(col) > max_gap_fraction:
            continue
        counter = Counter(col)
        counter.pop('-', None)
        if len(counter) < 2:
            continue
        most_common = counter.most_common(2)
        consensus += most_common[1][0]
    with open(output_file, 'w') as f:
        f.write(">next_best_consensus\n")
        f.write(consensus + "\n")
    

def main(blast_xml, pdb_file, output_dir="results"):
    """
    Run the transporter analysis pipeline.

    Parameters
    ----------
    blast_xml : str or Path
        BLASTP XML input file.
    pdb_file : str or Path
        PDB structure used to extract the query protein sequence.
    output_dir : str or Path
        Directory in which all pipeline outputs will be written.
    """
    blast_xml = Path(blast_xml).expanduser()
    pdb_file = Path(pdb_file).expanduser()
    output_dir = Path(output_dir).expanduser()

    if not blast_xml.is_file():
        raise FileNotFoundError(f"BLAST XML input file not found: {blast_xml}")
    if not pdb_file.is_file():
        raise FileNotFoundError(f"PDB input file not found: {pdb_file}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Retrieve protein sequences from all BLASTP hits.
    ret_fasta = output_dir / "retrieved.fasta"
    extract_and_fetch_sequences(blast_xml, ret_fasta)

    # Remove short sequences and sequences containing unknown residues.
    filtered_fasta = output_dir / "filtered.fasta"
    filter_sequences(ret_fasta, filtered_fasta)

    # Extract the query sequence from the supplied PDB structure.
    query_pdb_fasta = output_dir / "query_from_pdb.fasta"
    extract_sequence_from_pdb(pdb_file, query_pdb_fasta)

    # Filter retrieved sequences by coverage of the PDB-derived query.
    alignment_filtered = output_dir / "alignment_filtered.fasta"
    filter_by_alignment_coverage(
        query_pdb_fasta,
        filtered_fasta,
        alignment_filtered
    )

    # Cluster sequences and extract CD-HIT representative sequences.
    clustered_prefix = output_dir / "clustered"
    cluster_reps_fasta = output_dir / "cluster_reps.fasta"
    run_cdhit(alignment_filtered, clustered_prefix)
    parse_clusters(
        str(clustered_prefix) + ".clstr",
        alignment_filtered,
        cluster_reps_fasta
    )

    # Align representative sequences, trim the alignment, and generate consensus.
    aligned_fasta = output_dir / "aligned.fasta"
    cleaned_aln = output_dir / "cleaned_alignment.fasta"
    run_mafft(cluster_reps_fasta, aligned_fasta)
    filter_alignment(aligned_fasta, cleaned_aln)

    consensus_file = output_dir / "next_best_consensus.fasta"
    next_best_consensus(cleaned_aln, consensus_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the Ascomycota transporter protein analysis pipeline."
    )
    parser.add_argument(
        "blast_xml",
        help="Path to the BLASTP XML input file."
    )
    parser.add_argument(
        "--pdb",
        required=True,
        help="Path to the query protein PDB structure."
    )
    parser.add_argument(
        "--output-dir",
        default="results",
        help="Directory for pipeline outputs (default: results)."
    )

    args = parser.parse_args()
    main(args.blast_xml, args.pdb, args.output_dir)

