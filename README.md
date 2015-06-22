This repository contains the following:

- The Chicago R package (source and build)
- The PCHiCdata R data package with small example datasets for mouse and human 
- chicagoTools: scripts for preparing input files, running Chicago and processing the output  

To install the Chicago package, data package and chicagoTools (and all the R packages they depend on), download the whole repository and run:

```
#!bash

    Rscript setupChicago.R

```

In some cases, you may need to run the setup script with custom parameters:

```
#!bash

    Rscript setupChicago.R [--chicago-path=<Chicago-package-path>] [--data-path=<PCHiCdata-package-path>] [--rlib=<r-lib-dir>] [--bin=<chicagoTools-target-dir>]
    
```

The description of these parameters (all of them optional) is as follows:

 - chicago-path: the location of the Chicago package (either as a directory or as a tar.gz file). Defaults to the current directory
 - data-path: the location of the PCHiCdata package (either as a directory or as a tar.gz file). Defaults to the current directory
 - rlib: the location of the R library directory (known to R and for which you have write permissions). If not provided, the default R library directory is used (available via .libPaths()[1])
 - bin: if provided, chicagoTools will be moved to this path, alternatively they will be left at their current location. 

Note that the R packages require R version >= 3.1.2, and chicagoTools require bedtools, perl and python >= 2.7 that need to be pre-installed and added to PATH.

Alternatively, you can use standard R tools to only install the Chicago and/or PCHiCdata R packages. 

Please refer to the Chicago package vignette for more information. 

Chicago is developed and maintained by:

- Jonathan Cairns 
- Paula Freire Pritchett
- Steven Wingett
- Mikhail Spivakov ([spivakov@babraham.ac.uk](mailto:spivakov@babraham.ac.uk))

We are based at the [Regulatory Genomics Group](http://www.regulatorygenomicsgroup.org), [Babraham Institute](http://www.babraham.ac.uk), Cambridge UK.