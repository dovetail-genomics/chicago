**This is the repo for the second version of CHiCAGO capture-HiC peak calling pipeline** 

Currently CHiCAGO is still under development, so below are just some very rough instructions on how to install and run the pipeline in its current form.

Preparing restriction enzyme-specific data files (if you are not using HindIII)
=======================

**If you are using HindIII**, all the required files are available as a tarball archive at the following location:
https://www.dropbox.com/s/4nws1b3yqkfvwjd/HindIII_files.tar.gz?dl=0 (it's 221Mb so probably a bit too large for bitbucket)

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

Installing CHiCAGO
===============

1. Check that you have a recent version of bedtools (ours is v2.17.0) installed and executable as $ bedtools.

2. Open process_chic_single_core.sh and locate the ## EDIT ME line close to the top of the file. Hardcode the full path to the folder, in which *process_hicup.sh* is located, into the ‘pipelinedir=’ command immediately below this line. 

3. Open production_line_CHiCAGO2.R and edit the following:
 - Locate the ####EDIT ME#### line at the top of the file. Change the path to chicago.R and Functions_new_datatable.R, as will as the path to the resource file directory (fileDir variable).
 - Hard-code the full paths of the restriction enzyme-specific resource files under ### Resource file locations.
 - Check that parameter definitions under ### Fragment filtering and other settings match those used for generating restriction enzyme-specific resource files.

4. R dependencies.
 - Make sure you're running a recent version of R (ours is v3.0.3) 
 - CHiCAGO core requires the following R packages:

MASS
data.table
matrixStats
Hmisc
Delaporte

If some of them aren't currently installed, do so via install.packages()

  - The feature enrichment script requires Bioconductor and the add-on package GenomicRanges.
If not currently installed, do so: 
 source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges") 




Note for developers: Eventual aim is for the structure to be as: http://nvie.com/posts/a-successful-git-branching-model/
Currently, we are using the master branch & feature branches only. Upon the first release, a development branch will be added.

Developers: Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett.