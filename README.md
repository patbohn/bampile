# TLDR [WIP, not yet ready]

This is a small tool that helps quantifying linked mutations from sequencing data. 
It takes a bed-like file as input, where one can specify positions of interest, 
together with wildtype and (several) mutant alleles. The tool then iterates through
a bam file to count for each read the number of matches to the wildtype or mutant 
base, thus resulting in a match score. This is written out into a csv file. 

## Position of interest file format

The format to specify positions of interest is derived from the bed format, except
that there can be multiple columns for different mismatch patterns. 

## Installation

The software is a single binary that can be either downloaded from the Releases page, 
or alternatively (if it doesn't run on your system) you can compile it yourself. 

### Compilation 

This software is written in Rust and can thus be compiled with `cargo`. Follow the 
install instructions specific to your OS, then clone this github, `cd` into it and
run `cargo build --release`. The compiled binary will then be at `target/release/bampile`. 
