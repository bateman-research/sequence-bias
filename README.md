# Sequence composition and bias analyses

This repository contains scripts to analyse the sequence composition and bias of protein and nucleotide (DNA) sequences.

### Software requirements

Scripts are implemented in `R` v3.6, and use the following packages: `seqinr`, `argparse`, `ggplot2`, `dplyr` and `tidyr`.

### Scripts

This repository contains five scripts, listed below. 
Each script can be run from the command line with custom input and output options. 
Use the help option (-h/--help) on each script to see all the available options.

- [bias_profile_dna.R](bias_profile_dna.R): plot the nucleotide composition and bias profile along the DNA sequence of a gene.
- [bias_profile_prot.R](bias_profile_prot.R): plot the amino acid composition and bias profile along a protein sequence.
- [composition_dna.R](composition_dna.R): analyse the DNA composition for a subset of genes, including GC content, nucleotide asymmetries and expected vs observed nucleotide frequencies, and generate PCA and dendogram plots.
- [composition_prot.R](composition_prot.R): analyse the amino acid composition for a subset of protein sequences, and generate PCA and dendogram plots.
- [dotplot.R](dotplot.R): generate a dot-plot of two protein or DNA sequences.

### Results

Each script takes as input a single or multiple protein sequences in a FASTA format, and other options and parameters, and produces figures and/or results tables.
Results and figures for example input files are provided in the [examples](examples) folder.

### Reference

These scripts were created and used for my PhD research work at the European Bioinformatics Institute (EMBL-EBI).
More details are provided in my Thesis (insert link when available).

Aleix Lafita - April 2021
