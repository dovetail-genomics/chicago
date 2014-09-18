**This is the repo for the second version of CHiCAGO capture-HiC peak calling pipeline** 

Currently CHiCAGO is still under development, so below are just some very rough instructions on how to install and run the pipeline in its current form.

**0. Prepare restriction enzyme-specific data files (if you are not using HindIII)**

If you are using HindIII, all the required files are available as a tarball archive at the following location:
https://www.dropbox.com/s/4nws1b3yqkfvwjd/HindIII_files.tar.gz?dl=0 (it's 221Mb so probably a bit too large for bitbucket)

You will need to generate the following files yourself (the corresponding files available for HindIII are in brackets):

- A restriction fragment map (Digest_Human_HindIII.bed)
- Bait fragment map, using the same fragment IDs as the above map (Digest_Human_HindIII_baits.bed)
- Same with an additional column listing feature IDs for output (Digest_Human_HindIII_baits_ID.bed)

*Additionally, python scripts are provided to generate the following files:*

- Number of other ends per distance bin within the proximal range (Digest_Human_HindIII_NperBin.txt)

$ python countNperBin.py [--minFragLen=<min-restriction-fragment-size>] [--maxFragLen=<max-restriction-fragment-size] [--maxLBrownEst=<max-distance-for-estimating-brownian-noise>] 
                         [--binsize=<bin-size-in-bps-for-Brownian-noise-parameter-estimation>] [--removeAdjacent==<t|f>] 
						 --rmapfile=<restriction-fragment-map-filename> --baitmapfile=<bait-fragment-map-filename> --outfile=<output-filename> --picklefile=<output-python-pickle-file-name>

*- Number of baits per distance bin within the proximal range from the point of view of the other ends (Digest_Human_HindIII_NbaitsPerBin.txt)*

$ python countNbaitsPerBin.py [--maxLBrownEst=<max-distance-for-estimating-brownian-noise>] [--binsize=<bin-size-in-bps-for-Brownian-noise-parameter-estimation>] [--removeAdjacent==<t|f>] 
						 --rmapfile=<restriction-fragment-map-filename> --baitmapfile=<bait-fragment-map-filename> --outfile=<output-filename> --picklefile=<output-python-pickle-file-name>

*- List of all other ends falling into the proximal range for at least one bait (prox_OE.txt)*

$ python getProxOE.py [--minFragLen=<min-restriction-fragment-size>] [--maxFragLen=<max-restriction-fragment-size] [--maxLBrownEst=<max-distance-for-estimating-brownian-noise>] 
                      [--binsize=<bin-size-in-bps-for-Brownian-noise-parameter-estimation>] [--removeAdjacent==<t|f>] 
					  --rmapfile=<restriction-fragment-map-filename> --baitmapfile=<bait-fragment-map-filename> --outfile=<output-filename> --picklefile=<output-python-pickle-file-name>


**Important.** The input parameters for the above three scripts must match the corresponding parameters (with the same name) for the R CHiCAGO pipeline hard-coded into production_line_CHiCAGOv2.R. 
CHiCAGO will check whether they matched though and terminate with an error if not. The defaults for all optional parameters match those currently hard-coded into the production line script.





Note for developers: Eventual aim is for the structure to be as: http://nvie.com/posts/a-successful-git-branching-model/
Currently, we are using the master branch & feature branches only. Upon the first release, a development branch will be added.

Developers: Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett.

