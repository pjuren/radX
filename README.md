RADIX: Regression Analysis of Differential Inclusion of eXons 
=============================================================

RADIX is software for identifying exons which show differential 
rates of inclusion within a gene between different sample types. 
It supports multiple replicates per sample type, and arbitrary 
experimental design.  


Contact Information
-------------------

Philip J. Uren
uren@usc.edu
http://smithlabresearch.org/

Installation
------------
*Before attempting to compile RADX please make sure that GNU Scientific 
Library (http://www.gnu.org/software/gsl/) is installed on your system*
Alternatively, you can download pre-compiled binaries for either Lunux or Mac 
from http://smithlabresearch.org/software/

To compile RADX, enter the program's root directory (e.g. radmeth/) and  
execute

> make

After the compilation, the binaries can be found in radmeth/bin/

Usage
-----

radx takes as input a 


Input format
------------
Readcounts for exons are provided to RADIX in BED format files. 
The experimental design is provided in a tab-seperated text file.
The counts for each sample are provided in a BED file. This file
must contain the following tab-seperated columns: chromosome, exon
start, exon end, exon name, score (the count of reads in this exon),
and exon strand (+ or -). Here are some example entries:

chr1	1334022	1334345	exon1_exon_gene1	34	+
chr1  1334808 1334899 exon2_exon_gene1	10	+
chr3	124904	1250023	exon1_exon_gene2	24	-

exons must follow this naming style: exonName_exon_geneName, where
exonName and geneName can be any string desired, but must form
a unique identifier. All entries with the same value for 'geneName' 
are considered to be part of the same gene. The exons for a given 
gene must be non-overlapping. Read counts for each sample should be 
provided in seperate BED files. 

Output format
-------------
RADIX produces output in BED format. A single BED file is produced, 
regardless of how many samples are provided. 

License
-------
Copyright (C) 2013 University of Southern California and
							 Philip J Uren
               Egor Dolzhenko
               Andrew D Smith

    Authors: Philip J Uren, Andrew D. Smith and Egor Dolzhenko

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

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

