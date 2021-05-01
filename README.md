# sequence-bias

This repository contains scripts to analyse the sequence composition and bias of protein and nucleotide (DNA) sequences.

### Software requirements

Scripts are implemented in `R` v3.6, and use the following packages: `seqinr`, `argparse`, `dplyr` and `tidyr`.

### Scripts

This repository contains four scripts, listed below. 
Each script can be run from the command line with custom input and output options. 
Use the help option (-h/--help) on each script to see all the available options.

- [bias_profile_dna.R](bias_profile_dna.R): 
- [bias_profile_prot.R](bias_profile_prot.R): 
- [composition_dna.R](composition_dna.R): 
- [composition_prot.R](composition_prot.R): 
- [dotplot.R](dotplot.R)

### Results

Each script takes as input a single or multiple protein sequences in a FASTA file and other options, and produces several figures and results files.
Results and figures for example input files are provided in the [examples](examples) folder.

### Reference

These scripts were created and used for my PhD research work at the European Bioinformatics Institute (EMBL-EBI).
More details are provided in my Thesis (insert link when available).

Aleix Lafita - April 2021
