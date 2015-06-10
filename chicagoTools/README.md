chicagoTools are an assorted set of scripts associated with the CHiCAGO R package.  

Currently it includes the following software:

- Scripts for preparing "design files" needed for the CHiCAGO package:
    makeNBaitsPerBinFile.py
    makeNPerBinFile.py
    makeProxOEFile.py
    
- Script to processing BAM files into CHiCAGO input files:
    bam2chicago.sh
    
- The wrapper script for the CHiCAGO package itself:
    runChicago.R
    
- Scripts for post-processing CHiCAGO output:
    makePeakMatrix.R - used for bundling the calls from multiple samples into a single data matrix
    [FIXME] FitDistCurve.Rmd - used for reestimating CHiCAGO p-walue weighting parameters based on user data
    
