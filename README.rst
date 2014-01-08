Introduction: what is fragrecruit?
---------------------------------

**fragrecruit** is a pure-Python application for generating fragment
recruitment plots from metagenomic data.  It operates on the output from
fragment recruitment tools such as FR-HIT or any tool that provides output
in the SAM format (such as Bowtie or BWA) and generates plots that describe
how sequencing reads align to reference genomes.

Prerequisites
-------------

**fragrecruit** has been tested with Python version 2.6.  matplotlib must
also be installed.

Installing
----------

The source distribution for the most recent version can be obtained from
the `fragrecruit project page <https://github.com/mstrosaker/fragrecruit>`_.
All of the python files in the source directory are needed, as well as
a configuration file; see the sample.config file to see how this config file
should be formatted.

How to use it?
--------------

Run "./fragplot -h" for a usage message:

Usage::
  ./fragplot [options] <config_file>

  Options
    -v: verbose
    -h: display this help message
    -o <filename>: string to prepend the output file name
    -m <n>: only plot sequence alignments of length n or greater
    -i <n>: minimum identity to include in plot (y-axis extent)
    -g <taxonids>: only generate plots for the specified genomes
       multiple genomes can be specified in a comma-separated list
    -s <samples>: only plot points for the specified samples
       multiple samples can be specified in a comma-separated list
    -b <n>: begin the plot at the nth nucleotide in the ref genome
    -l <n>: include n nucleotides in the x-axis of the plot
       the default behavior is to plot the entirety of the genome

Specify the genomic data for the reference genomes and the alignment data
in the file passed as <config_file>; see sample.config for an example.

License
-------

Copyright (c) 2014 Michael Strosaker.  See the LICENSE file for license
rights and limitations (MIT).

