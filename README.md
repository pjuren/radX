RADIX: Regression Analysis of Differential Inclusion of eXons 
=============================================================

RADIX is software for identifying exons which show differential 
rates of inclusion within a gene between different sample types. 
It's core model is built around a beta-binomial regression, allowing
it to support multiple replicates per sample type, and arbitrary 
experimental design. 

Copyright and License Information
---------------------------------
Copyright (C) 2014
University of Southern California,
Philip J. Uren, Andrew D. Smith, Egor Dolzhenko 
  
Authors: Philip J. Uren, Andrew D. Smith, Egor Dolzhenko
    
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
      
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
        
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

This software package contains Google Test and BAMTools -- see their 
respective directories for copyright and license information for those 
packages. 

Building and Installing RADIX
-----------------------------
RADIX has been designed to operate in a UNIX-like environment.
It has been tested on MacOS X Mavericks, USCLinux, Ubuntu and SUSE. We've
made every effort to ensure it's compatible with a wide range of systems, but
we can't test on all possible platforms and configurations. If you find an
incompatibility when building or using it, please let us know!

### Step 0 -- dependencies and requirements ###

  64-bit machine and GCC version > 4.1 (to support TR1). 
    
  RADIX additionally requires a functioning installation of the GNU 
  Scientific Library (GSL). If you don't already have this installed, you will
  need to download and install it from http://www.gnu.org/software/gsl/
  This release has been tested with GSL version 1.14
  Your installation of GSL must be built for a 64 bit architecture.
                
  [OPTIONAL] Tests
    The regression tests (which you don't need to run, but can to test that 
    RADIX was built correctly) require Python. If you don't have it, you can 
    download Python from http://www.python.org/
    This release has been tested with Python 2.6.1 
                                    
  
### Step 1 -- configuring ###

  First configure the installation. To do this, where '>' is your prompt and 
  the CWD is the root of the distribution, type:
      
  > ./configure 

### Step 2 -- building ###
  
  To build the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution  
        
  > make all 
            
### Step 3 -- installing ###
  
  To install the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution
                  
  > make install
                      
  This will place the binaries in the bin directory under the package root.
  They can be used directly from there without any additional steps. You can
  add that directory to your PATH environment variable to avoid having to 
  specify their full paths, or you can copy the binaries to another directory
  of your choice in your PATH.
                                  
### Step 4 -- testing [OPTIONAL] ###
  
  You can verify that RADIX was built correctly by running the included
  regression tests. At the command prompt (assuming your prompt is '>') 
  type the following:
                                          
  > make test 

Usage
-----

### Basic usage ###

RADIX takes as input a set of read counts for exons in each sample
and a design matrix describing the experimental setup. These two inputs
are always required, and in the simplest use of the program nothing else
need be specified. For example, if your design matrix was in the file 
`design.mat` and you had four replicates (say two for some experimental 
condition and two control replicates) with read counts in 
`control_rep_1.bed`, `control_rep_2.bed`, `exp_rep_1.bed` and 
`exp_rep_2.bed` then you would run the program as follows

> radix -d design.mat control_rep_1.bed control_rep_2.bed exp_rep_1.bed exp_rep_2.bed

Note that if you have more than one condition, you must also specify the 
option '-f <condition>' to indicate which condition you wish to test. 
See below for more details about conditions and how to specify them in 
the design matrix. Additional options can be specified to adjust, for 
example the signficance threshold used by RADIX. Run the program with 
no arguments to see a list of other options that can be used:

> radix

### Exon read-count file format ###

Readcounts for exons are provided to RADIX in BED format files. 
The counts for each sample are provided in a seperate file. These files
must contain the following tab-seperated columns: chromosome, exon
start, exon end, exon name, score (the count of reads in this exon),
and exon strand (+ or -). Here are some example entries:

```
chr1  1334022 1334345 gene1_exon_1  34  +
chr1  1334808 1334899 gene1_exon_2  10  +
chr3  124904  1250023 gene2_exon_1  24  -
```

exons must follow this naming style: geneName_exon_exonName, where
exonName and geneName can be any string desired, but must form
a unique identifier. All entries with the same value for 'geneName' 
are considered to be part of the same gene. The exons for a given 
gene must be non-overlapping, and all input files must contain 
read-counts for the same set of exons. If either of those conditions
is not met, the program will exit with an error. BED files are tab
seperated. 

### Design matrix format ###

The first row of the design matrix must contain each of the experimental
conditions in your data, expressed in such a way that you can answer the 
question of whether that condition is present or not in each sample. As
a simple example, consider testing two control replicates against two 
cancer replicates. Here, there is just one condition: cancer, and each
sample either has, or does not have this condition (note that there is 
no need to include "control" as a condition since it is just the inverse
of the "cancer" condition). Each subsequent line in the file should 
start with a filename that contains the exon read counts for a sample
and be followed by either a zero or a one indicating whether the 
corresponding condition is present for that sample or not (i.e. the first 
column after the name corresponds to the first condition listed in the 
file's first line, the second column corresponds to the second condition 
listed in the file's first line, etc.) 
Here's an example:

```
                     some_experimental_condition
control_rep_1.bed                0
control_rep_2.bed                0
exp_rep_1.bed                    1
exp_rep_2.bed                    1
```

Exact alignment and additional whitespace does not need to be maintained,
we use it here just for display purposes. 

### Output format ###

RADIX produces output in BED format. A single BED file is produced, 
regardless of how many samples are provided. Each exon present in the input
files will also be present in the output. The score field of the BED 
file (the fifth column) will be the mean log fold-change for the exon's 
inclusion rate. There will be an addition field (beyond the default six
columns in a BED file) that gives the p-value the resulted from testing 
whether the exon was differentially included or not. Only exons with a 
p-value less than the significance threshold (set by default to 0.01, 
but adjustable) are output.

Contacts and Bug Reports
------------------------
Philip J. Uren
uren@usc.edu

Egor Dolzhenko 
dolzhenk@usc.edu

Andrew D. Smith
andrewds@usc.edu

If you found a bug in RADIX, we'd like to know about it. Before contacting us
though, please check the following list:

1.) Are you using the latest version? The bug you found may already have 
    been fixed.
2.) Check that your input is in the correct format and you have selected
    the correct options.
3.) Please reduce your input to the smallest possible size that still 
    produces the bug; we will need your input data to reproduce the 
    problem.

