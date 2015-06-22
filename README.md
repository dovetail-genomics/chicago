**CHiCAGO: Capture HiC Analysis of Genomic Organisation** 

CHiCAGO is a set of tools for Capture HiC data analysis. A preprint describing the statistical algorithm behind Chicago interaction calling will be available in the next few weeks. 

This repository contains the following files:

- The Chicago R package     
- The PCHiCdata R data package with small example Promoter Capture HiC datasets for mouse and human  
- chicagoTools: scripts for preparing input files, running Chicago and processing the output  

Please refer to the Chicago R package vignette and the chicagoTools README file for more information.

Standard R tools (such as ```R CMD INSTALL``` and ```install.packages()```) can be used to install the Chicago and PCHiCdata packages. The chicagoTools scripts require no specific installation.   
   
Note that the R packages require R version >= 3.1.2, and chicagoTools require bedtools, perl and python >= 2.7 that need to be pre-installed and added to PATH.   
   
Alternatively, to install the Chicago package, PCHiC data package and chicagoTools (as well as all R packages they depend on), download the whole repository and run:

```
#!bash

    Rscript setupChicago.R

```

In some cases, you may need to run setupChicago.R with custom parameters:

```
#!bash

    Rscript setupChicago.R [--chicago-path=<Chicago-package-path>] [--data-path=<PCHiCdata-package-path>] [--rlib=<r-lib-dir>] [--bin=<chicagoTools-target-dir>]
    
```

The description of these parameters (all of them optional) is as follows:

 - chicago-path: the location of the Chicago package (either as a directory or as a tar.gz file). Defaults to the current directory
 - data-path: the location of the PCHiCdata package (either as a directory or as a tar.gz file). Defaults to the current directory
 - rlib: the location of the R library directory (known to R and for which you have write permissions). If not provided, the default R library directory is used (available via .libPaths()[1])
 - bin: if provided, chicagoTools will be moved to this path, alternatively they will be left at their current location. 

Chicago is developed and maintained by:

- Jonathan Cairns 
- Paula Freire Pritchett
- Steven Wingett
- Mikhail Spivakov ([spivakov@babraham.ac.uk](mailto:spivakov@babraham.ac.uk))

We are based at the [Regulatory Genomics Group](http://www.regulatorygenomicsgroup.org), [Babraham Institute](http://www.babraham.ac.uk), Cambridge UK.