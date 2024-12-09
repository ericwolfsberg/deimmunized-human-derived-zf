# Code for "Machine-Guided Dual-Objective Protein Engineering for Deimmunization and Therapeutic Functions"
This repository contains code for generating human-derived deimmunized zinc finger arrays, in zf_arrays_for_site.py
and jupyter notebooks for replicating code-derived results that are included in figures 5b, 5c, 5d, and 6b.

## ZF Targeting
The intended usage for zf_arrays_for_site.py is
```
python3 zf_arrays_for_site.py -t <DNA target sequence from 5' to 3'> -o <output filename> -b <binding affinity data> -m <MHC-acceptable pair data>  [-s <1 to use ZF scores, 0 not to; default 1>, -n <number of arrays to print, or 0 to print all; default 0>]
```
The -b flag can take values of Z or D for precalculated ZifRC and DeepZF affinity data for human ZFs, or a filename for custom user-generated data. 
Custom affinity data files should be tab-separated with columns containing the name of each ZF, its amino acid sequence, its target codon, and the score if scoring is to be used.
Any rows with fewer than four tab-separated columns if the -s flag is 1 or three columns if it is 0 will be ignored.


The -m flag can take values of M or N for precalculated MARIA and NetMHCIIpan acceptable ZF pairs, or a filename for custom user-generated data. 
Custom ZF pair files should contain all acceptable ZF-ZF pairs separated by tabs.
Any ZF names which are not present in the affinity data file will be ignored.


If target sequence is not a multiple of three nucleotides it will be truncated at the 3' end
### Requirements
* python > 3.6

## Data reproduction notebooks
The Jupyter notebooks within the jupyter and data replication contain code for replicating the data within the indicated figures. For the NetMHCII site analyses in Fig 5d Data.ipynb please unzip the files within the "netmhc final pairs" directory.
