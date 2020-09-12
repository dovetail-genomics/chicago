Tuning CHiCAGO settings for four-cutter enzymes

Valeriya Malysheva, Helen Ray-Jones and Mikhail Spivakov

The following default settings of the CHiCAGO pipeline have been optimised for Promoter Capture Hi-C experiments using the 6-cutter enzyme HindIII:

minFragLen = 150,
maxFragLen = 40000,
maxLBrownEst = 1.5e6,
binsize=20000,
weightAlpha = 34.1157346557331, 
weightBeta = -2.58688050486759,
weightGamma = -17.1347845819659,
weightDelta = -7.07609245521541

When running CHiCAGO with four-cutter enzymes (or other frequent cutters), these parameters should be adjusted. As a starting point we recommend using the following parameters:

minFragLen = 75,
maxFragLen = 1200,
maxLBrownEst = 75000,
binsize=1500,
weightAlpha = 24.5, 
weightBeta = -2.16,
weightGamma = -21.2,
weightDelta = -9.2

However, the optimal values of these parameters depend on factors beyond the choice of the restriction enzyme and in particular, on sequencing depth. It is therefore advisable to assess the suitability of specific settings for the analysed data, using plotBackgroundSparsity.R in chicagoTools.
The four-cutter weight parameters are conservative and we recommend that users estimate optimal weights for their experimental settings using fitDistCurve.R in chicagoTools.

When adjusting the parameters, it is important to bear in mind that CHiCAGO scores from interaction calls generated with different parameter settings may not be directly comparable.
