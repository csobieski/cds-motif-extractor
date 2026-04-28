# sequence-extractor

Python tool to extract genomic features from GFF3 files based on DNA or protein motif matches.

This repository contains a Python script for extracting annotated genomic features (e.g., CDS, rRNA, tRNA, ncRNA, genes) from GFF3 files containing embedded FASTA sequences. The tool searches for user-defined DNA or protein motifs (including regular expressions) and outputs matching regions in FASTA format. It is particularly suited for bacterial genome analyses but can be applied to any GFF3-formatted dataset.

## Features
- Extract any GFF3 feature type (CDS, rRNA, tRNA, gene, ncRNA, etc.)
- Search motifs in DNA or translated protein sequences
- Support for regular expression searches
- Handles multi-exon / joined features
- Optional upstream padding
- Compatible with bacterial translation tables

## Requirements
- Python 3
- biopython
- gffutils

Install dependencies:
```bash
pip install -r requirements.txt

## Example usage using query of DNA sequence:

python extract_cds_by_motif.py \
  --input-dir "Genome_KP" \
  --query "ATGCGT" \
  --output extracted_features.fasta \
  --search-in dna \
  --feature-type CDS \
  --padding 400

## Example usage using query of protein sequence:

python extract_cds_by_motif.py \
  --input-dir "Genome_KP" \
  --query "MKTIIALSY" \
  --output extracted_protein_hits.fasta \
  --search-in protein \
  --feature-type CDS
