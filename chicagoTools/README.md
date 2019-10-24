chicagoTools are an assorted set of scripts associated with the Chicago R package.  

Currently, the following software is included:

- The script for preparing "design files" needed for the Chicago package:
       + makeDesignFiles.py

- Deprecated individual scripts for preparing the same "design files":  
       + makeNBaitsPerBinFile.py  
       + makeNPerBinFile.py  
       + makeProxOEFile.py  
          
- The script for processing BAM files into Chicago input files:  
       + bam2chicago.sh  
   
- The wrapper script for the Chicago package:  
       + runChicago.R  
    
- The script for bundling the interaction calls from multiple samples into a single data matrix:  
       + makePeakMatrix.R  
       
- The script for reestimating Chicago p-value weighting parameters based on user data:
       + fitDistCurve.R  
  
**The script for preparing the "design files" needed for the Chicago package**
  
A new set of design files needs to be generated for each combination of restriction enzyme and captured baits (i.e., when
.baitmap and/or .rmap files are updated). (Note that it is possible to also create "virtual digests" in the same way, 
where the actual restriction fragments are pooled into larger "virtual" fragments). The python script makeDesignFiles.py prepares all of these file. The script os compatible with python version 2.7+ (but not 3.x).  

**Important**: The makeDesignFiles.py script needs to be re-run when updating the experimental design with respect to either restriction enzyme, capture design or Chicago settings with respect to binning, proximal distance range and/or restriction fragment filtering. 
  
The script uses as input the following files:  
  
- rmap file (```.rmap```): a tab-separated file of the format ```<chr> <start> <end> <numeric ID>```, describing the restriction digest (or "virtual digest" 
  if pooled fragments are used). These numeric IDs are referred to as "otherEndID" in Chicago. All fragments mapping outside of the digest coordinates will be disregarded 
  by both these scripts and Chicago.   
- baitmap file (```.baitmap```): a tab-separated file of the format ```<chr> <start> <end> <numeric ID> <annotation>```, listing the coordinates of the 
  baited/captured restriction fragments (should be a subset of the fragments listed in rmapfile), their numeric IDs (should match those listed in 
  rmapfile for the corresponding fragments) and their annotations (such as, for example, the names of baited promoters). The numeric IDs are referred to as "baitID" in Chicago.    
  
The files it generates are ASCII files containing the following information:  
  
- NPerBin file (.npb): <baitID> <Total no. valid restriction fragments in distance bin 1> ... <Total no. valid restriction fragments in distance bin N>,  
where the bins map within the "proximal" distance range from each bait (0; maxLBrownEst] and bin size is defined by the binsize parameter.  
- NBaitsPerBin file (.nbpb): <otherEndID> <Total no. valid baits in distance bin 1> ... <Total no. valid baits in distance bin N>,   
where the bins map within the "proximal" distance range from each other end (0; maxLBrownEst] and bin size is defined by the binsize parameter.  
- Proximal Other End (ProxOE) file (.poe): <baitID> <otherEndID> <absolute distance>     
for all combinations of baits and other ends that map within the "proximal" distance range from each other (0; maxLBrownEst].  
Data in each file is preceded by a comment line listing the input parameters used to generate them.  
   
The script takes the following input parameters:  
   
```python makeDesignFiles.py [--designDir=.] [--rmapfile=designDir/*.rmap] [--baitmapfile=designDir/*.baitmap]  [--outfilePrefix=designDir/<rmapfileName>] --minFragLen=150 --maxFragLen=40000 [--maxLBrownEst=1500000] [--binsize=20000] [--removeb2b=True] [--removeAdjacent=True]```
   
- The following parameters specify the input files that need to be created prior to running these scripts and the output file name:  
    + rmapfile: path to rmap file   
    + baitmapfile: path to baitmap file        
    + outfilePrefix: the name of the output file prefix, including the path (if not provided, will have the same name as the rmap file and the extension specific to each of the
  three file types: .npb, .nbpb and .poe, respectively.   
    + designDir: if rmapfile, baitmapfile and/or outfile are not explicitly specified, the script will automatically look for rmapfile and baitmapfile at this location (under the extensions 
.rmap and .baitmap, respectively), and will also place the output files there.   
          
- The following options should be consistent with the corresponding settings in the Chicago R package, and the script needs to be rerun whenever these settings are modified:   
    + binsize: the size of the bins (in bps) for pooling restriction fragments   
    + minFragLen: the min fragment length cutoff (no default: recommended 75 for DpnII, 150 for HindIII)
    + maxFragLen: the max fragment length cutoff (no default: recommended 1200 for DpnII, 40000 for HindIII)
    + maxLBrownEst: the "proximal distance range" for estimating Brownian noise    
   
- The following options should always be set to defaults in the current implementation of Chicago:  
    + removeb2b: True, meaning that bait-to-bait interactions should not be counted when computing the total numbers of fragments at a given distance.  
    + removeAdjacent: True, meaning that fragments immediately adjacent to bait should not be counted.  

Note that minFragLen and max no longer have defaults to avoid common mistakes.

**The script for processing BAM files into Chicago input files**

The Unix shell script bam2chicago.sh takes as input a BAM file corresponding to aligned Capture HiC reads for a single experiment and the capture design files .rmap and .baitmap.  
The output is an ASCII .chinput file used as input by the Chicago package. 

The BAM file is a paired-end file produced by a HiC aligner; Chicago has only been tested with data produced by HiCUP (http://www.bioinformatics.babraham.ac.uk/projects/hicup/). However, it should theoretically be possible to use 
other HiC aligners for this purpose. 

The formats of the design files (.baitmap and .rmap) are as described above:  
- rmap file (```.rmap```): a tab-separated file of the format ```<chr> <start> <end> <numeric ID>```, describing the restriction digest (or "virtual digest" 
  if pooled fragments are used). These numeric IDs are referred to as "otherEndID" in Chicago. All fragments mapping outside of the digest coordinates will be disregarded 
  by both these scripts and Chicago.   
- baitmap file (```.baitmap```): a tab-separated file of the format ```<chr> <start> <end> <numeric ID> <annotation>```, listing the coordinates of the 
  baited/captured restriction fragments (should be a subset of the fragments listed in rmapfile), their numeric ids (should match those listed in 
  rmapfile for the corresponding fragments) and their annotations (such as, for example, the names of baited promoters). The numeric IDs are referred to as "baitID" in Chicago.    

The output ```.chinput``` file has the format ```<baitID> <otherEndID> <N> <otherEndLen> <distSign>```, where N is the number of reads detected for ligation products between the "bait" and "other end", otherEndLen 
is the length of the "other-end" restriction fragment and distSign is the linear distance between the bait and other-end fragments, respectively.   

In addition, bam2chicago.sh produces a paired-end bed (```.bedpe```) file that lists all read pairs corresponding to bait-to-bait interactions. This file is not used for Chicago; the rationale for creating it is that 
bait-to-bait interactions are symmetric and may also be analysed by any HiC normalisation/interaction calling tool (in addition to Chicago) if desired.
    
The shell script requires a bash-compatible Unix shell, perl and bedtools. Bedtools can be obtained at https://github.com/arq5x/bedtools2; they need to be installed and added to ```$PATH```.  
  
The script takes the following input parameters:   
   
```bam2chicago.sh <bamfile> <baitmap-file> <digest-rmap-file> <sample-name> [nodelete]```

- bamfile: path to the input BAM file  
- baitmap-file: path to the input .baitmap file
- digest-map-file: path to the input .rmap file
- sample-name: will be used for naming the output folder, in which the output files will be placed and as the basename of the output .chinput and .bedpe files  
- nodelete: a flag to prevent the script from deleting intermediate files  
  
**The wrapper script for the Chicago package**   
   
The R script runChicago.R can be used to run a typical Chicago analysis. Please refer to Chicago vignette and the inline help for the individual R functions used
for more details on each analysis step.   
   
runChicago.R performs the following steps:   
   
- Creates the chicagoData object given the design folder and, if needed, with other custom settings using ```setExperiment()```   
- Reads in the input file(s) and merges replicates if necessary using ```readSample()``` or ```readAndMerge()```    
- Runs interaction calling using ```chicagoPipeline()```
- Saves the full chicagoData object as an R image (in the Rds or RDa format)   
- Exports significant interactions in a genome browser-readable format using ```exportResults()```
- Plots the profiles multiple random baits using ```plotBaits()```   
- Estimates the enrichment of significant interactions for user-specified genomic 
features versus distance-matched controls using ```peakEnrichment4Features()```   
- Saves the settings used for running Chicago and the input parameters of the script itself to a text file   
- Sorts output files into a directory tree with the subfolders ```data/```, ```diag_plots/```,  ```enrichment_data/``` and ```examples/```   
   
A typical analysis will use the following options:   
   
```Rscript runChicago.R --design-dir DESIGN-DIR --en-feat-list EN-FEAT-LIST <input-files> <output-prefix>```   
   
- design-dir: the name of the design folder containing exactly one of the following file types: 
```.baitmap```, ```.rmap```, ```.npb```, ```.nbpb``` and ```.poe``` with the corresponding extensions (see above and the 
Chicago package documentation for the description of these formats). The option defaults to the current directory.   
- en-feat-list: the name of the feature list file of the format <feature-name> <feature-bed-file-location>, where feature-bed-files  
contain the genomic coordinates of the features to compute the enrichment at Chicago significant interactions 
(Chicago signal cutoff 5 is used by default).   
- <input-files> - a single .chinput file (produced by ```bam2chicago.sh``` or any other method) or a comma-separated 
list of .chinput files corresponding to the multiple biological replicates of the same experimental condition. Note that technical replicates 
should instead be deduplicated and pooled prior to running bam2chicago and submitted as a single .chinput file.
- <output-prefix> - experiment name used in the naming of the output folder and as a prefix for output file names.
   
The full list of the parameters is however more extensive:    
   
```Rscript runChicago.R [--help] [--print-memory] [--rda] [--save-df-only] [--examples-full-range] [--settings-file SETTINGS-FILE] 
[--design-dir DESIGN-DIR] [--cutoff CUTOFF] [--export-format EXPORT-FORMAT] [--export-order EXPORT-ORDER] [--examples-prox-dist EXAMPLES-PROX-DIST] [--output-dir OUTPUT-DIR] 
[--en-feat-list EN-FEAT-LIST] [--en-feat-files EN-FEAT-FILES] [--en-feat-folder EN-FEAT-FOLDER] [--en-full-cis-range] [--en-trans] 
[--en-min-dist EN-MIN-DIST] [--en-max-dist EN-MAX-DIST] [--en-sample-no EN-SAMPLE-NO] [--features-only]
<input-files> <output-prefix>
```   
   
- help: print a help message
- print-memory: print memory use during ```chicagoPipeline()``` execution   
- rda: save the image of the chicagoData object as an RDa (under the name cd) rather than an Rds file   
- save-df-only: save the image of only the data table part of the chicagoData object (```cd@x```) as a data frame (i.e., converting from 
data.table to data frame and discarding the ```@params``` and ```@settings``` slots)
- examples-full-range: in addition to plotting interactions within 1Mb from baits, also plot the same for the full distance range   
- settings-file: the path to a settings file, from which to load custom ```chicagoData@settings```   
- design-dir: the name of the design folder containing exactly one of the following file types: 
```.baitmap```, ```.rmap```, ```.npb```, ```.nbpb``` and ```.poe``` with the corresponding extensions (see above and the 
Chicago package documentation for the description of these formats). The option defaults to the current directory.   
- cutoff: a signal cutoff to use for significant interactions [default: 5]   
- export-format: file format for writing out peaks: one or more of the following: seqMonk,interBed,washU_text,washU_track (comma-separated) 
[default: washU_text]   
- export-order: should the results be ordered by "score" or genomic "position"? [default: position]   
- examples-prox-dist: the distance limit for plotting "proximal" examples [default: 1Mb]   
- en-feat-list: the name of the feature list file of the format <feature-name> <feature-bed-file-location>, where feature-bed-files  
contain the genomic coordinates of the features to compute the enrichment at Chicago significant interactions 
(Chicago signal cutoff 5 is used by default).   
- en-feat-files: a comma-separated list of files with genomic feature coordinates for computing peaks' enrichment (to provide them explicitly 
instead of using the --en-feat-list option)   
- en-feat-folder: the folder, in which all feature files are located (if provided, --en-feature-file(s) don't need to list the full path)   
- en-full-cis-range: assess the enrichment for features for the full distance range (same chromosome only; use --en-trans in addition to include trans-interactions). Can be very slow!   
- en-trans: include trans-interactions into enrichment analysis   
- en-min-dist: the lower distance limit for computing enrichment for features (default: 0)   
- en-max-dist: the upper distance limit for computing enrichment for features (default: 1Mb)   
- en-sample-no: the number of negative samples, over which to compute the enrichment for features (default: 100)   
- features-only: re-run feature enrichment analysis with Chicago output files. 
With this option, <input-files> must be either a single Rds file (must contain full Chicago objects) or '-', 
in which case the file location will be inferred automatically from <output-prefix> and files added to the corresponding folder.   
- <input-files> - a single .chinput file (produced by bam2chicago.sh or any other method) or a comma-separated list of file names 
corresponding to the multiple biological replicates of the same experimental condition. Note that technical replicates should 
instead be deduplicated and pooled prior to running bam2chicago and submitted as a single .chinput file.   
- <output-prefix> - experiment name used in the naming of the output folder and as a prefix for output file names.   
      
**The script for bundling the interaction calls from multiple samples into a single data matrix**  
   
When running multiple samples through CHiCAGO it is convenient to represent the results in the form of a "peak matrix". This matrix lists the coordinates, annotations and sample-wise scores for all interactions that pass a signal threshold in at at least one sample. The peak matrix can then be used for downstream analyses such as clustering by interaction and sample type and integration with other types of data.   
   
The R script makePeakMatrix.R takes as input the list of chicago output data images (by default, the Rds files containing the chicagoData objects) and outputs a peak matrix as a text and Rds file. In addition the script generates a hierarchical clustering dendrogram of the samples based on the peak matrix scores.   
   
A typical run of ```makePeakMatrix.R``` will use the following options:  
    
```Rscript makePeakMatrix.R [--twopass] <names-file> <output-prefix>```   
   
- <names-file>: full path to a tab-separated file with sample names (1st column) and full paths to input Rds files (2nd column)   
- <output-prefix>: the prefix to use for the output files   
- twopass: first obtain a list of significant interactions as a union of peaks in each dataset, then reload and subset for these interactions prior to merging. Slower but uses significantly less memory.     
   
The full list of options for ```makePeakMatrix.R``` is listed below:   
   
```Rscript makePeakMatrix.R [--help] [--twopass] [--notrans] [--vanilla] [--rda] [--var VAR] [--print-memory] [--scorecol SCORECOL] [--cutoff CUTOFF] [--fetchcol FETCHCOL] [--lessthan] [--maxdist MAXDIST] [--digestmap DIGESTMAP] [--baitmap BAITMAP] [--peaklist PEAKLIST] [--clustmethod CLUSTMETHOD] [--clustsubset CLUSTSUBSET] <names-file> <output-prefix>```    
   
- help: print a help message   
- twopass: first obtain a list of significant interactions as a union of peaks in each dataset, then reload and subset for these interactions prior to merging. Slower but uses significantly less memory.   
- notrans: exclude trans-chromosomal interactions from the peak matrix (an alternative, but faster way to save memory compared with ```--twopass```)   
- vanilla: a flag indicating that the input RDa/RDS images contain only the  data frames and not Chicago objects. In this case, the ```--digestmap``` and ```--baitmap``` options are required (see below).   
- rda: load data from an RDa archive rather than the default Rds. In this case, the name of the variable in the RDa containing the image is given by the ```--var``` option (see below).   
- var: the name of the variable containing the ```chicagoData``` object or the peak data frame in the RDa images (default: x)   
- print-memory: print memory info at each step   
- scorecol: the column name in the ```chicagoData@x``` slot or the input data frame containing the Chicago scores (default: score). Note that this also allows to create peak matrices for entities other than Chicago scores (e.g., the raw or normalised reads).   
- cutoff: the ```scorecol``` signal cutoff   
- fetchcol: Instead of collecting scores for the peak matrix, choose the name of a different column to collect information from. ```scorecol``` is still used to threshold interactions. Note that if ```fetchcol``` is different from ```scorecol```, the ```--twopass``` mode will be enforced. 
- lessthan: pick interactions with ```scorecol``` below the cutoff rather than above    
- maxdist: max distance from bait to include into the peak matrix   
- digestmap: full path to digest map file; will override settings from ```chicagoData``` if provided. Required for ```--vanilla```.   
- baitmap: full path to bait map ID file; will override settings from ```chicagoData``` even if provided. Required for ```--vanilla```.   
- peaklist: use a predefined peak list (such as the one generated by the first pass of the ```--twopass``` mode)   
- clustmethod: the clustering method to use for clustering columns (average/ward.D2/complete) (default: average)   
- clustsubset: number of interactions to randomly subset for clustering. Full dataset used if total number of interactions in the peak matrix is below this number. (default: 1e+06)   
   
**The script for reestimating Chicago p-value weighting parameters based on user data**  

CHiCAGO uses a p-value weighting procedure to upweight proximal interactions and downweight distal interactions. To weight appropriately, CHiCAGO needs to know how the probability of an interaction event between two fragments decreases as the distance between these fragments increases. The parametrization used has four parameters: alpha, beta, gamma, and delta (more details on the exact parametrization are given in the CHiCAGO paper).  

We have calibrated these parameters on high-confidence calls from seven human Macrophage data sets (i.e. interactions that pass our p-value threshold in all seven samples). Provided that your cell type is not too dissimilar to these calibration data, it should be fine to leave the parameters at their default settings. However, if your data set is from an unusual cell type, you may wish to recalibrate these parameters using data from cell types similar to yours. This script provides a CHiCAGO .settings file, which contains estimates of each of the four parameters, and a plot to show how the interaction abun.  

Users must specify <output-prefix>, and either inputs or summaryInput. Typically, the first time you run this script, you specify --inputs:  

Rscript fitDistCurve.R cellType --inputs 1stFile.Rda,2ndFile.Rda,3rdFile.Rda  

This procedure generates a file called cellType_summaryInput.Rda. To rerun the same inputs but with different parameters, you can save time by specifying --summaryInput instead:  

Rscript fitDistCurve.R cellTypeLargerBin --summaryInput cellType_summaryInput.Rda --largeBinSize 2000000  

```Rscript fitDistCurve.R [--help] [--opts OPTS] [--inputs INPUTS] [--summaryInput SUMMARYINPUT] [--threshold THRESHOLD] [--subsets SUBSETS] [--largeBinSize LARGEBINSIZE] [--binNumber BINNUMBER] [--halfNumber HALFNUMBER] <output-prefix>```  

- <output-prefix>: All output files will begin with the value of this argument.  
- inputs: Comma-separated list, specifying locations of saved chicagoData objects, in either .Rda or .Rds form. If .Rda, there must be only one chicagoData object in the file.  
- summaryInput: An .Rda file of summary information -- the max P-val for each putative interaction, and the location of the .rmap file. This file will be generated if it wasn't provided.  
- threshold: Threshold applied to log(p) values (NB: not the CHiCAGO score!). In other words, log(p) must be below this value in all of the samples. [default: -10]  
- subsets: To ensure robustness, the data are partitioned into approximately equal subsets. Parameters are estimated separately on each subset. Then, an estimate for each parameter is derived from its median value across subsets. The number of subsets is controlled by this argument. [default: 5]  
- largeBinSize, binNumber, halfNumber: Parameters pertaining to the bins used in the analysis. Default breaks occur at 0, 31.25k, 62.5k, 125k, 250k, 500k, 1m, 2m, 3m, 4m, ..., 16m. The breaks are constructed by taking [binNumber] bins of size [largeBinSize], then breaking the first bin into two, iterating [halfNumber] times. [default: 1000000, 16, 5]
