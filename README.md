**CHiCAGO: Capture HiC Analysis of Genomic Organisation** 

CHiCAGO is a set of tools for Capture HiC data analysis. 

CHiCAGO is presented in [this paper](http://www.genomebiology.com/2016/17/1/127): Cairns J / Freire-Pritchett P, Wingett SW, V rnai C, Dimond A, Plagnol V, ZerbinoÿD, Schoenfelder S, Javierre B-M, Osborne C, Fraser P, Spivakov M. CHiCAGO: Robust Detection of DNA Looping Interactions in Capture Hi-C data. Genome Biology. 2016. 17:127. 

This repository contains the following files:

- The Chicago R package     
- The PCHiCdata R data package with small example Promoter Capture HiC datasets for mouse and human  
- chicagoTools: scripts for preparing input files, running Chicago and processing the output  

Please refer to the Chicago R package vignette and the chicagoTools README file for more information.

*Compatibility notices*

- CHiCAGO is not compatible with R package Delaporte v2.3.0. The compatibility issue has been kindly fixed but Delaporte package author in all subsequent versions.

- CHiCAGO is currently not compatible with bedtools v2.26 due to BED format compliance checking introduced in this version. Please do not upgrade from v2.25 while we are working to resolve this issue.

*News*

- A new version of bam2chicago script is released, bam2chicago_V02.sh, with improved usability and compatibility with HiCUP combinations script. 

- Suggested parameter set for four-cutters (tested with DpnII): maxLBrownEst=75000, binsize=1500, minFragLen=75, maxFragLen=1200.

- Check out Chicdiff, our new differential caller for Capture Hi-C data that works jointly with Chicago. Chicdiff is available on [github](https://github.com/RegulatoryGenomicsGroup/chicdiff/) and is presented in [this paper](https://doi.org/10.1093/bioinformatics/btz450): Cairns J / Orchard W / Malysheva V, Spivakov M. Chicdiff: a computational pipeline for detecting differential chromosomal interactions in Capture Hi-C data. Bioinformatics. 2019. AOP: btz450. 

- Version 1.13: Default values of parameters related to restriction fragment sizes and binning are now propagated automatically from the design files. Therefore they only need to be provided to makeDesignFiles.py and not the Chicago package separately (as with previous versions).

- chicagoTools: makeDesignFiles.py no longer has defaults for minFragLen and maxFragLen to avoid mistakes. Our recommended settings are: for HindIII - 150 and 40000, respectively; for DpnII - 75 and 1200, respectively. 

- Version 1.1.5: Default values of tlb.minProxOEPerBin and tlb.minProxB2BPerBin have changed. No action is required unless you specified non-default values, or wish to re-run the pipeline on old chicagoData objects. See the [NEWS](https://bitbucket.org/chicagoTeam/chicago/src/master/Chicago/NEWS?fileviewer=file-view-default) file for more details.

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

The R packages are also part of Bioconductor starting from version 3.3, and installation using ```biocLite()``` is available. However, as Bioconductor releases only happen twice a year, more recent versions of the R packages may be available from here.

If you encounter any problems, please [post an issue](https://bitbucket.org/chicagoTeam/chicago/issues) or email the developers. In the email, include output from the R command ``sessionInfo()``, along with any error messages encountered.

*Contact information*

Chicago was originally developed by Jonathan Cairns, Paula Freire Pritchett, Steven Wingett and Mikhail Spivakov ([mikhail.spivakov@lms.mrc.ac.uk](mailto:mikhail.spivakov@lms.mrc.ac.uk)), with recent modifications and parameter tuning by Valeriya Malysheva, Helen Ray-Jones and Monica Della Rosa.


CHiCAGO was developed at the Regulatory Genomics Group, Babraham Institute, Cambridge UK. From July 2018, the group is based at MRC London Institute of Medical Sciences in London, where it is known as [Functional Gene Control](http://www.lms.mrc.ac.uk/groups/functional-gene-control) group.

More details (including the full credits) can be found at [www.functionalgenecontrol.group/chicago](http://www.functionalgenecontrol.group/chicago).
