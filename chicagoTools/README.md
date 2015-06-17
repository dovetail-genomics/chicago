chicagoTools are an assorted set of scripts associated with the Chicago R package.  

Currently it includes the following software:

- Scripts for preparing "design files" needed for the Chicago package:  
       + makeNBaitsPerBinFile.py  
       + makeNPerBinFile.py  
       + makeProxOEFile.py  
    
- The script for processing BAM files into Chicago input files:  
       + bam2chicago.sh  
   
- The wrapper script for the Chicago package:  
       + runChicago.R  
    
- The script for bundling the interaction calls from multiple samples into a single data matrix:  
       + makePeakMatrix.R  
       
- The script for reestimating Chicago p-walue weighting parameters based on user data:
       + [FIXME] FitDistCurve.Rmd  
  
**Scripts for preparing the "design files" needed for the Chicago package**
  
A new set of design files needs to be generated for each combination of restriction enzyme and captured baits (i.e., when
.baitmap and/or .rmap files are updated). (Note that it is possible to also create "virtual digests" in the same way, 
where the actual restriction fragments are pooled into larger "virtual" fragments). Each python script (makeNPerBinFile.py, makeNBaitsPerBinFile.py,
makeProxOEFile.py) prepares one type of design file, respectively. They are compatible with python version 2.7+ (but not 3.x).  

**Important**: All of these scripts need to be re-run (in any order or in parallel) when updating the experimental design with respect to  
either restriction enzyme, capture design or Chicago settings with respect to binning, proximal distance range and/or restriction fragment filtering. 
  
These scripts take the following input files:  
  
- rmap file (.rmap): a tab-separated file of the format <chr> <start> <end> <numeric ID>, describing the restriction digest (or "virtual digest" 
  if pooled fragments are used). These numeric IDs are referred to as "otherEndID" in Chicago. All fragments mapping outside of the digest coordinates will be disregarded 
  by both these scripts and Chicago.   
- baitmap file (.baitmap): a tab-separated file of the format <chr> <start> <end> <numeric ID> <annotation>, listing the coordinates of the 
  baited/captured restriction fragments (should be a subset of the fragments listed in rmapfile), their numeric ids (should match those listed in 
  rmapfile for the corresponding fragments) and their annotations (such as, for example, the names of baited promoters). The numeric IDs are referred to as "baitID" in Chicago.    
  
These files are ASCII files containing the following information:  
  
- NPerBin file (.npb): <baitID> <Total no. valid restriction fragments in distance bin 1> ... <Total no. valid restriction fragments in distance bin N>,  
where the bins map within the "proximal" distance range from each bait (0; maxLBrownEst] and bin size is defined by the binsize parameter.  
- NBaitsPerBin file (.nbpb): <otherEndID> <Total no. valid baits in distance bin 1> ... <Total no. valid baits in distance bin N>,   
where the bins map within the "proximal" distance range from each other end (0; maxLBrownEst] and bin size is defined by the binsize parameter.  
- Proximal Other End (ProxOE) file (.poe): <baitID> <otherEndID> <absolute distance>     
for all combinations of baits and other ends that map within the "proximal" distance range from each other (0; maxLBrownEst].  
Data in each file is preceded by a comment line listing the input parameters used to generate them.  
   
All three scripts take the same input parameters:  
   
``` python makeNPerBinFile.py / makeNBaitsPerBinFile.py / makeProxOEFile.py  
	[--designDir=.]  
    [--rmapfile=designDir/*.rmap]  
	[--baitmapfile=designDir/*.baitmap]  
	[--outfile=designDir/<rmapfileName>.<npb|nbpb|poe>]  
    [--minFragLen=150]   
    [--maxFragLen=40000]   
    [--maxLBrownEst=1500000]  
    [--binsize=20000]  
    [--removeb2b=True]  
    [--removeAdjacent=True] ```
   
- The following parameters specify the input files that need to be created prior to running these scripts and the output file name:  
    + rmapfile: path to rmap file   
    + baitmapfile: path to baitmap file        
    + outfile: the name of the output file (if not provided, will have the same name as the rmap file and the extension specific to each of the
  three file types: .npb (makeNPerBinFile.py), .nbpb (makeNBaitsPerBinFile.py) and .poe (makeProxOEFile.py), respectively.   
    + designDir: if rmapfile, baitmapfile and/or outfile are not explicitly specified, the scripts will automatically look for rmapfile and baitmapfile at this location (under the extensions 
.rmap and .baitmap, respectively), and will also place the output files there.   
          
- The following options should be consistent with the corresponding settings in the Chicago R package, and the scripts need to be rerun whenever these settings are modified:   
    + binsize: the size of the bins (in bps) for pooling restriction fragments   
    + minFragLen: the min fragment length cutoff    
    + maxFragLen: the max fragment length cutoff   
    + maxLBrownEst: the "proximal distance range" for estimating Brownian noise    
   
- The following options should always be set to defaults in the current implementation of Chicago:  
    + removeb2b: True, meaning that bait-to-bait interactions should not be counted when computing the total numbers of fragments at a given distance.  
    + removeAdjacent: True, meaning that fragments immediately adjacent to bait should not be counted.  


** The script for processing BAM files into Chicago input files **

The Unix shell script bam2chicago.sh takes as input a BAM file corresponding to aligned Capture HiC reads for a single experiment and the capture design files .rmap and .baitmap.  
The output is an ASCII .chinput file used as input by the Chicago package. 

The BAM file is a paired-end file produced by a HiC aligner; Chicago has only been tested with data produced by HiCUP (http://www.bioinformatics.babraham.ac.uk/projects/hicup/). However, it should theoretically be possible to use 
other HiC aligners for this purpose. 

The formats of the design files (.baitmap and .rmap) are as described above:
- rmap file (.rmap): a tab-separated file of the format <chr> <start> <end> <numeric ID>, describing the restriction digest (or "virtual digest" 
  if pooled fragments are used). These numeric IDs are referred to as "otherEndID" in Chicago. All fragments mapping outside of the digest coordinates will be disregarded 
  by both these scripts and Chicago.   
- baitmap file (.baitmap): a tab-separated file of the format <chr> <start> <end> <numeric ID> <annotation>, listing the coordinates of the 
  baited/captured restriction fragments (should be a subset of the fragments listed in rmapfile), their numeric ids (should match those listed in 
  rmapfile for the corresponding fragments) and their annotations (such as, for example, the names of baited promoters). The numeric IDs are referred to as "baitID" in Chicago.    

The output .chinput file has the format <baitID> <otherEndID> <N> <otherEndLen> <distSign>, where N is the number of reads detected for ligation products between the "bait" and "other end", otherEndLen 
is the length of the "other-end" restriction fragment and distSign is the linear distance between the bait and other-end fragments, respectively.   

In addition, bam2chicago.sh produces a paired-end bed (.bedpe) file that lists all read pairs corresponding to bait-to-bait interactions. This file is not used for Chicago; the rationale for creating it is that 
bait-to-bait interactions are symmetric and may also be analysed by any HiC normalisation/interaction calling tool (in addition to Chicago) if desired.
    
The shell script requires a bash-compatible Unix shell, perl and bedtools. Bedtools can be obtained at https://github.com/arq5x/bedtools2; they need to be installed and added to $PATH.  
  
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
file corresponding to the multiple biological replicates of the same experimental condition.
Note that technical replicates should instead be deduplicated and pooled prior to running bam2chicago and submitted as a single .chinput file.
- <output-prefix> - experiment name used in the naming of the output folder and as a prefix for output file names.
   
The full list of the  parameters is however much more extensive:    
   
```Rscript runChicago.R [--help] [--print-memory] [--rda] [--save-df-only] [--examples-full-range] [--settings-file SETTINGS-FILE] 
[--design-dir DESIGN-DIR] [--cutoff CUTOFF] [--export-format EXPORT-FORMAT] [--export-order EXPORT-ORDER] [--examples-prox-dist EXAMPLES-PROX-DIST] [--output-dir OUTPUT-DIR] 
[--en-feat-list EN-FEAT-LIST] [--en-feat-files EN-FEAT-FILES] [--en-feat-folder EN-FEAT-FOLDER] [--en-full-cis-range] [--en-trans] 
[--en-min-dist EN-MIN-DIST] [--en-max-dist EN-MAX-DIST] [--en-sample-no EN-SAMPLE-NO] [--features-only]
<input-files> <output-prefix>```   
   
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
