This repository contains the following:

- The Chicago R package (source and build)
- The PCHiCdata R data package with small example datasets for mouse and human 
- chicagoTools: scripts for preparing input files, running Chicago and processing the output  

To install the Chicago package, data package and chicagoTools (and all the R packages they depend on), download the whole repository and run:

```
#!bash

    Rscript setupChicago.R [--bin=<scripts-target-dir>] [--path=<chicago-package-path>] [--rlib=<r-lib-dir>] 

```

Note that the R packages require R version >= 3.1.2, and chicagoTools require bedtools, perl and python >= 2.7 that need to be pre-installed and added to PATH.

Alternatively, you can use standard R tools to only install the Chicago and/or PCHiCdata R packages. 

Please refer to the Chicago package vignette for more information. 

Chicago is developed and maintained by Jonathan Cairns, Paula Freire Pritchett, Steven Wingett and Mikhail Spivakov (mailto:spivakov@babraham.ac.uk) at the Regulatory Genomics Group, Babraham Institute, Cambridge UK (http://www.regulatorygenomicsgroup.org ; http://www.babraham.ac.uk).