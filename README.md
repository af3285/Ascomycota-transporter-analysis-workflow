# Ascomycota Transporter Analysis Workflow

## Overview

This repository contains a Python-based bioinformatics workflow for identifying and analyzing candidate transporter proteins in Ascomycota. It combines sequence retrieval, sequence filtering, clustering, multiple sequence alignment, and consensus sequence generation.

The workflow uses BLASTP to identify candidate protein sequences, retrieves corresponding sequences from NCBI, and compares them against a reference protein sequence extracted from a supplied PDB structure. Candidate sequences are filtered by length, residue composition, and alignment coverage before being clustered with CD-HIT. Representative sequences are then aligned using MAFFT, and the resulting alignment is processed to generate a next-best consensus sequence.

A Bash wrapper (`Pipeline_runner`) is also included to automate the initial BLASTP search and execute the Python analysis pipeline.

This repository showcases the computational workflow and its implementation. Original research datasets and input files are not included.

## Workflow

The analysis consists of the following steps:

1. **BLASTP search:** Identify candidate protein homologs using BLASTP. The wrapper searches the NCBI `nr` database remotely and restricts results to Ascomycota.
2. **Sequence retrieval:** Parse BLASTP XML output, extract unique protein accession identifiers, and retrieve corresponding sequences from NCBI using Entrez.
3. **Initial sequence filtering:** Remove sequences shorter than 300 amino acids or containing unknown amino acid residues (`X`).
4. **Reference extraction:** Extract the longest polypeptide sequence from the supplied PDB structure.
5. **Coverage filtering:** Perform local pairwise alignments against the reference sequence and retain candidates covering at least 90% of the reference sequence.
6. **Sequence clustering:** Cluster candidate sequences using CD-HIT at 80% sequence identity and 90% shorter-sequence coverage, retaining cluster representatives.
7. **Multiple sequence alignment:** Align representative sequences using MAFFT with the L-INS-i strategy.
8. **Alignment filtering:** Remove alignment columns containing more than 20% gaps.
9. **Consensus generation:** Generate a next-best consensus sequence by selecting the second-most frequent amino acid at eligible alignment positions. Positions with fewer than two distinct non-gap residues are skipped.

## Repository Structure

| File | Description |
|---|---|
| `AnalysisCode_configurable.py` | Main Python analysis pipeline |
| `Pipeline_runner` | Bash wrapper for BLASTP search and pipeline execution |
| `README.md` | Project documentation |

## Requirements

### Software

- Python 3
- [Biopython](https://biopython.org/)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
- [CD-HIT](https://github.com/weizhongli/cdhit)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)

The wrapper uses remote BLASTP, so a locally installed copy of the NCBI `nr` database is not required. An internet connection is required for remote BLASTP searches and NCBI Entrez sequence retrieval.

### Python dependency installation

Install Biopython using pip:

```bash
python3 -m pip install biopython
```

BLAST+, CD-HIT, and MAFFT must also be installed separately and accessible through the system PATH.

### NCBI Entrez configuration

The Python pipeline uses Biopython's Entrez interface to retrieve protein sequences from NCBI. The `Entrez.email` setting in `AnalysisCode_configurable.py` should be updated to an appropriate email address before execution.

## Running the Workflow

There are two ways to execute the analysis.

### Option 1: Run the complete workflow using the Bash wrapper

The wrapper automates the BLASTP search and subsequent Python analysis.

The wrapper expects:
- A protein FASTA file supplied as a command-line argument.
- A reference PDB file named `6uwf.pdb` in the working directory.

From the repository directory, make the script executable:

```bash
chmod +x Pipeline_runner
```

Then execute it with a FASTA input:

```bash
./Pipeline_runner query.fasta
```

The wrapper performs the following operations:

1. Submits a remote BLASTP search against the NCBI `nr` database.
2. Restricts the search to Ascomycota using the `Ascomycota[Organism]` Entrez query.
3. Applies an E-value threshold of `1e-5`.
4. Saves the results in XML format as `blast_output.xml`.
5. Executes `AnalysisCode_configurable.py` using the BLASTP results and `6uwf.pdb`.
6. Writes the analysis outputs to the `results/` directory.

**Important:** The wrapper currently uses fixed filenames for the BLAST output and reference PDB. It also assumes it is executed from the repository directory. Users should ensure that `6uwf.pdb` is present before running it.

### Option 2: Run the Python pipeline using existing BLASTP results

The Python script can also be executed independently when a BLASTP XML file has already been generated.

The pipeline requires two input files:

- A BLASTP XML file containing protein sequence hits, generated using `-outfmt 5`.
- A PDB structure file containing the reference protein sequence.

Execute the pipeline with:

```bash
python3 AnalysisCode_configurable.py blast_output.xml --pdb reference.pdb --output-dir results
```

The output directory is optional and defaults to `results/`. It is created automatically if it does not already exist.

## Output Files

The pipeline generates intermediate files and a final consensus sequence in the specified output directory.

| File | Description |
|---|---|
| `retrieved.fasta` | Protein sequences retrieved from NCBI |
| `filtered.fasta` | Sequences passing the initial length and residue filters |
| `query_from_pdb.fasta` | Reference protein sequence extracted from the PDB structure |
| `alignment_filtered.fasta` | Sequences passing the reference coverage filter |
| `clustered` | CD-HIT clustered sequence output |
| `clustered.clstr` | CD-HIT cluster membership information |
| `cluster_reps.fasta` | Representative sequences selected from CD-HIT clusters |
| `aligned.fasta` | Multiple sequence alignment generated by MAFFT |
| `cleaned_alignment.fasta` | Alignment after gap-based column filtering |
| `next_best_consensus.fasta` | Generated next-best consensus sequence |

## Configuration

Several analysis parameters can be adjusted directly in `AnalysisCode_configurable.py`:

| Parameter | Default | Description |
|---|---:|---|
| Minimum sequence length | 300 amino acids | Minimum length for retained protein sequences |
| Reference coverage | 90% | Minimum fraction of the reference sequence covered by a candidate |
| CD-HIT identity | 80% | Sequence identity threshold for clustering |
| CD-HIT coverage | 90% | Minimum coverage of the shorter sequence during clustering |
| Maximum alignment gap fraction | 20% | Maximum fraction of sequences allowed to contain a gap at an alignment position |

The wrapper's BLASTP E-value threshold is set to `1e-5`.

## Input Data and Reproducibility

The original research input files are not included in this repository. Users must provide their own protein FASTA query and compatible PDB structure. Alternatively, users can provide an existing BLASTP XML file when running the Python pipeline independently.

The workflow retrieves protein sequences directly from NCBI, meaning that results may vary as database contents and annotations change over time.

The pipeline is intended as a research and exploratory bioinformatics workflow. Its output should be interpreted in the context of the biological question and the quality of the input sequences and reference structure.
