**CHiCAGO: Capture HiC Analysis of Genomic Organisation** 

CHiCAGO is a set of tools for Capture HiC data analysis. A preprint describing the statistical algorithm behind CHiCAGO is available on bioRxiv: [http://dx.doi.org/10.1101/028068](http://dx.doi.org/10.1101/028068).

This repository contains the following files:

- The Chicago R package     
- The PCHiCdata R data package with small example Promoter Capture HiC datasets for mouse and human  
- chicagoTools: scripts for preparing input files, running Chicago and processing the output  

Please refer to the Chicago R package vignette and the chicagoTools README file for more information.

*Installation instructions*

1. Make sure that you have R version >= 3.1.2. chicagoTools requires some additional dependencies: bedtools, perl and python >= 2.7 need to be pre-installed and added to PATH, plus the R package ``argparser`` is required - install with the R code:

```{r}
install.packages("argparser")
```

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
(This strategy downloads the repository multiple times. To avoid this, you can manually install Chicago's dependencies, then install the packages from the source directories using ``R CMD INSTALL`` or ``install.packages()``.)

The R packages have also been submitted to Bioconductor, and installation using ```biocLite()``` will be available upon acceptance.

If you encounter any problems, please [post an issue](https://bitbucket.org/chicagoTeam/chicago/issues) or email the developers. In the email, include output from the R command ``sessionInfo()``, along with any error messages encountered.

*Contact information*

Chicago is mainly developed and maintained by:

- Jonathan Cairns 
- Paula Freire Pritchett
- Steven Wingett
- Mikhail Spivakov ([spivakov@babraham.ac.uk](mailto:spivakov@babraham.ac.uk))

We are based at the [Regulatory Genomics Group](http://www.regulatorygenomicsgroup.org), [Babraham Institute](http://www.babraham.ac.uk), Cambridge UK.

More details (including the full credits) can be found at [regulatorygenomicsgroup.org/chicago](http://www.regulatorygenomicsgroup.org/chicago).