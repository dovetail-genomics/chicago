**CHiCAGO: Capture HiC Analysis of Genomic Organisation** 

CHiCAGO is a set of tools for Capture HiC data analysis. A preprint describing the statistical algorithm behind Chicago interaction calling will be available in the next few weeks. 

This repository contains the following files:

- The Chicago R package     
- The PCHiCdata R data package with small example Promoter Capture HiC datasets for mouse and human  
- chicagoTools: scripts for preparing input files, running Chicago and processing the output  

Please refer to the Chicago R package vignette and the chicagoTools README file for more information.

*Installation instructions*

1. Make sure that you have R version >= 3.1.2. In addition, chicagoTools require bedtools, perl and python >= 2.7 that need to be pre-installed and added to PATH.

2. Get the chicagoTools scripts by downloading the repository (downloads tab on left-hand menu). These scripts do not need to be installed further - see chicagoTools/README.md.

3. Install the R packages. An easy way to do this is by using functionality in devtools - run the following R code:
```{r}
install.packages("devtools")
library(devtools)
install_bitbucket("chicagoTeam/Chicago", subdir="Chicago")
```
Optionally, install the PCHiCdata package at the same time:
```{r}
install_bitbucket("chicagoTeam/Chicago", subdir="PCHiCdata")
```
(This strategy downloads the repository multiple times. To avoid this, you can manually install the packages from the source directories using ``R CMD INSTALL`` or ``install.packages()``.)

If you encounter any problems, please [post an issue](https://bitbucket.org/chicagoTeam/chicago/issues) or email the developers. In the email, include output from the R command ``sessionInfo()``, along with any error messages encountered.

*Contact information*

Chicago is developed and maintained by:

- Jonathan Cairns 
- Paula Freire Pritchett
- Steven Wingett
- Mikhail Spivakov ([spivakov@babraham.ac.uk](mailto:spivakov@babraham.ac.uk))

We are based at the [Regulatory Genomics Group](http://www.regulatorygenomicsgroup.org), [Babraham Institute](http://www.babraham.ac.uk), Cambridge UK.
