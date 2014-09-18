**This is the repo for the second version of CHiCAGO capture-HiC peak calling pipeline** 

Currently CHiCAGO is still under development, so below are just some very rough instructions on how to install and run the pipeline in its current form.

Preparing restriction enzyme-specific resource files
=======================

**If you are using HindIII**, all the required files are available as a tarball archive at the following location:
https://www.dropbox.com/s/l15wy6srherfuvn/HindIII_files.tar.gz?dl=0 (it's 221Mb so probably a bit too large for bitbucket)

**Otherwise**, you will need to generate the following files yourself (the corresponding files available for HindIII are in brackets):

- A restriction fragment map (Digest_Human_HindIII.bed)
- Bait fragment map, using the same fragment IDs as the above map (Digest_Human_HindIII_baits.bed)
- Same with an additional column listing feature IDs for output (Digest_Human_HindIII_baits_ID.bed)

**Additionally, if not using HindIII, python scripts are provided to generate the following files:**

- Number of other ends per distance bin within the proximal range (Digest_Human_HindIII_NperBin.txt)

$ python countNperBin.py [--minFragLen=<min-restriction-fragment-size>] [--maxFragLen=<max-restriction-fragment-size] [--maxLBrownEst=<max-distance-for-estimating-brownian-noise>] 
                         [--binsize=<bin-size-in-bps-for-Brownian-noise-parameter-estimation>] [--removeAdjacent==<t|f>] 
						 --rmapfile=<restriction-fragment-map-filename> --baitmapfile=<bait-fragment-map-filename> --outfile=<output-filename> --picklefile=<output-python-pickle-file-name>

- Number of baits per distance bin within the proximal range from the point of view of the other ends (Digest_Human_HindIII_NbaitsPerBin.txt)

$ python countNbaitsPerBin.py [--maxLBrownEst=<max-distance-for-estimating-brownian-noise>] [--binsize=<bin-size-in-bps-for-Brownian-noise-parameter-estimation>] [--removeAdjacent==<t|f>] 
						 --rmapfile=<restriction-fragment-map-filename> --baitmapfile=<bait-fragment-map-filename> --outfile=<output-filename> --picklefile=<output-python-pickle-file-name>

- List of all other ends falling into the proximal range for at least one bait (prox_OE.txt)

$ python getProxOE.py [--minFragLen=<min-restriction-fragment-size>] [--maxFragLen=<max-restriction-fragment-size] [--maxLBrownEst=<max-distance-for-estimating-brownian-noise>]  [--binsize=<bin-size-in-bps-for-Brownian-noise-parameter-estimation>] [--removeAdjacent==<t|f>] 
--rmapfile=<restriction-fragment-map-filename> --baitmapfile=<bait-fragment-map-filename> --outfile=<output-filename> --picklefile=<output-python-pickle-file-name>

**Important.** The input parameters for the above three scripts must match the corresponding parameters (with the same name) for the R CHiCAGO pipeline hard-coded into production_line_CHiCAGOv2.R. 

CHiCAGO will check whether they matched though and terminate with an error if not. The defaults for all optional parameters match those currently hard-coded into the production line script.

Installing CHiCAGOv2
===============

1. Check that you have a recent version of bedtools (ours is v2.17.0) installed and executable as $ bedtools.

2. Open process_chic_single_core.sh and locate the ## EDIT ME line close to the top of the file. Hardcode the full path to the folder, in which *process_hicup.sh* is located, into the ‘pipelinedir=’ command immediately below this line. 

3. Open run_chicago2.sh and locate the ##EDIT ME line at the top. Hardcode the full path to folder, in which *production_line_CHiCAGO2.R* is located, into the 'chicagopath=' command immediately below this line.

4. Open production_line_CHiCAGO2.R and edit the following:
 - Locate the ####EDIT ME#### line at the top of the file. Change the path to chicago.R and Functions_new_datatable.R, as will as the path to the resource file directory (fileDir variable).
 - Hard-code the full paths of the restriction enzyme-specific resource files under ### Resource file locations.
 - Check that parameter definitions under ### Fragment filtering and other settings match those used for generating restriction enzyme-specific resource files.

5. R dependencies.
 - Make sure you're running a recent version of R (ours is v3.0.3) 
 - CHiCAGO core requires the following R packages:

MASS, data.table, matrixStats, Hmisc, Delaporte

If some of these packages aren't currently installed, do so via install.packages()

  - The feature enrichment script requires Bioconductor and the add-on package GenomicRanges.
If not currently installed, do so: source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges") 

Running CHiCAGOv2
==========

**CHiCAGO is run in two steps:** 
-	process_chic_single_core.sh uses shell commands and bedtools to convert HiCUP BAM files into CHiCAGO input files.
-	production_line_v2.R  (that can be executed via the shell script run_fullchr_norep.sh) runs CHiCAGO and checks for enrichment of CHiCAGO’s significant ‘other ends’ (i.e., putative enhancers) for genomic regions of interest, such as, for example, ENCODE or BLUEPRINT regions.
The first step takes ~1.5h to run, the second one takes longer and requires ~20Gb RAM, so both are best run on a powerful machine (such as a cluster node).

*Instructions below are given assuming the standard HindIII resource files. Change the filenames accordingly if using another enzyme and custom files.*

1. From the folder where the HiCUP BAM file is located, run the shell script to convert BAM files to text files for R input:
----------
$ CHiCAGO/chic_tools/process_chic_single_core.sh <BAM file> <path>/Digest_Human_HindIII_baits.bed Digest_Human_HindIII.bed <name>
**Attention**: Digest_Human_HindIII_baits.bed, not Digest_Human_HindIII_baits_ID.bed file should be provided as parameter.
The output folder will be ./sample_<name>.
The file that that analysis pipeline will use is: 
    ./sample_<name>/<name>_bait_otherEnd_N_len_distSign.txt

2. Prepare a feature folder for your cell type (or use a related cell type).
-----------
For example, download from Blueprint ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/ several ChIP-seq datasets in bed format. 
In this folder, put a features.txt tab-delimited file of the format:
<dataset name>	<ChIP-seq bed file name>

3. From the folder where the HiCUP BAM file is located, run:
------------

**If have only one replicate:**

$ <path-to-CHiCAGOv2>/run_chicago2.sh sample_<name>/<name>_bait_otherEnd_N_len_distSign.txt <results-folder> <sample-name> <feature-folder-path> <full-path-to-features.txt file>

**If have two replicates:**

$ ./run_fullchr_norep.sh %2 <results-folder>/ <sample-name> < feature-folder-path>/ <full-path-to-features.txt file> sample_<name1>/<name1>_bait_otherEnd_N_len_distSign.txt sample_<name2>/<name2>_bait_otherEnd_N_len_distSign.txt

Note that you don't need to create any folders beforehand in either step 1 or step 3. 

**For <n> replicates**, use same syntax as above (but %<n> instead of %2), listing all sample files at the end of the command line.

CHiCAGO output
==========

The results will be in the ./<results-folder>/data, and various plots in the other subfolders of ./<results-folder>. 

**In the /data folder**, the .ibed and .txt files and are readable, respectively by WashU browser (epigenomegateway.wustl.edu) and Seqmonk (a 2-row format, where the first row corresponds to the bait and the second row to the respective other end). The score threshold of –log(adjusted p-value) of 11 is applied, but it is rather arbitrary.

**The data frame x** used to produce these files (one row per interaction) is stored in the .RDa file in the same folder.

The /examples folder contains PDFs with bait profiles for 25 random baits - for full chromosome length and zoomed in to 1Mb, respectively.

**The /overlap_plots folder** contains barplots showing the numbers of “enhancers” overlapping with genomic features of choice (yellow bars) versus the expected numbers computed using randomly sampled fragments chosen such that their distribution of distances from promoters matches that for the “enhancers”. Currently these overlaps are computed for interactions within 1 Mb from their respective baits only.

**In the /diag_plots folder**, we currently store diagnostic plots for other end normalisation and technical noise estimation, but more will be added.

Good luck!
----------
Just to repeat, this is a very early development version, so problems installing and running it are expected. Do not hesitate to contact Mikhail when this happens. 

Developers
-----------
Mikhail Spivakov (mikhail.spivakov@babraham.ac.uk), Jonathan Cairns, Paula Freire Pritchett.

Note for developers
----------
Eventual aim is for the structure to be as: http://nvie.com/posts/a-successful-git-branching-model/
Currently, we are using the master branch & feature branches only. Upon the first release, a development branch will be added.
