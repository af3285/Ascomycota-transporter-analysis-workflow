# Ascomycota Transporter Analysis Workflow

## Overview

This repository contains a Python-based bioinformatics workflow developed to analyze protein sequences associated with Ascomycota transporters. The pipeline integrates sequence retrieval, filtering, clustering, multiple sequence alignment, and consensus sequence generation to facilitate comparative protein sequence analysis.

The workflow uses BLASTP XML output to retrieve candidate protein sequences from NCBI, extracts a reference sequence from a supplied PDB structure, and applies sequence coverage filtering to identify relevant homologs. The resulting sequences are clustered using CD-HIT, aligned using MAFFT, and processed to generate a next-best consensus sequence.

This repository is intended as a showcase of the computational workflow and its implementation. The original input data and research files are not included.

## Workflow

1. **Sequence retrieval:** Extract protein accession identifiers from BLASTP XML results and retrieve corresponding sequences from NCBI.
2. **Sequence filtering:** Remove sequences below a specified length threshold or containing unknown amino acid residues.
3. **Coverage filtering:** Extract a reference protein sequence from a PDB structure and filter candidate sequences based on alignment coverage.
4. **Sequence clustering:** Use CD-HIT to group similar sequences and identify representative sequences.
5. **Multiple sequence alignment:** Align representative sequences using MAFFT.
6. **Alignment processing:** Remove alignment columns exceeding a specified gap threshold.
7. **Consensus generation:** Generate a next-best consensus sequence from the processed alignment.

## Requirements

- Python 3
- [Biopython](https://biopython.org/)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (for generating BLASTP XML input)
- [CD-HIT](https://github.com/weizhongli/cdhit)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)

The pipeline also requires an email address configured for NCBI Entrez access.

## Running the Pipeline

### 1. Prepare the input files

The pipeline requires two input files:

- A BLASTP XML output file containing protein sequence hits.
- A PDB structure file containing the reference protein sequence.

**Note:** These input files are not included in this repository. Users must supply their own data.

The BLASTP XML file should be generated using the XML output format (`-outfmt 5`). For example:

```bash
blastp -query query.fasta -db protein_database -out blast_output.xml -outfmt 5
```

This command assumes that a suitable BLAST protein database has already been installed or configured.

### 2. Configure NCBI access

Before running the pipeline, configure the email address used for NCBI Entrez requests in `AnalysisCode_configurable.py` by updating the `Entrez.email` setting.

### 3. Execute the pipeline

Run the script from the command line, specifying the BLASTP XML file, reference PDB structure, and desired output directory:

```bash
python3 AnalysisCode_configurable.py blast_output.xml --pdb reference.pdb --output-dir results
```

The output directory will be created automatically if it does not already exist.

### 4. Output files

The pipeline generates the following intermediate and final files:

| File | Description |
|---|---|
| `retrieved.fasta` | Protein sequences retrieved from NCBI |
| `filtered.fasta` | Sequences passing initial filtering |
| `query_from_pdb.fasta` | Reference sequence extracted from the PDB structure |
| `alignment_filtered.fasta` | Sequences passing query coverage filtering |
| `clustered` | CD-HIT clustered sequence output |
| `clustered.clstr` | CD-HIT cluster membership information |
| `cluster_reps.fasta` | Representative sequences selected from clusters |
| `aligned.fasta` | Multiple sequence alignment |
| `cleaned_alignment.fasta` | Alignment after gap-based filtering |
| `next_best_consensus.fasta` | Generated next-best consensus sequence |

## Notes

- The repository does not include the original BLASTP results, PDB structure, or other research input data.
- Users must provide their own compatible input files to execute the workflow.
- Filtering thresholds and analysis parameters can be adjusted within the Python script.
- The workflow is designed for protein FASTA sequences and assumes that the supplied BLASTP XML and PDB structure are compatible with the analysis steps.
