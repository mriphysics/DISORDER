# DISORDER
Tools for aligned reconstructions of volumetric MRI data

This repository provides tools to implement the methods and reproduce the experiments included in the manuscript ''Motion corrected MRI with DISORDER: Distributed and Incoherent Sample Orders for Reconstruction Deblurring using Encoding Redundancy'', L Cordero-Grande, G Ferrazzi, RPAG Teixeira, J O'Muircheartaigh, AN Price, and JV Hajnal, arXiv:1910.00540, 2019.

The code has been developed in MATLAB and has the following structure:

###### ./
contains the functions to run the illustrations and experiments included in Figs. 3-8 of the manuscript: *fig[0304,0506,0708].m*.

###### ./Acquisition
contains template functions to implement the DISORDER reorderings in MRI scanners: *electrostaticRepulsionDISORDER.m*, *samplingDISORDER.m*.

###### ./IO
contains functions for IO, visualization, and parameter setting: *disorderAlgorithm.m*, *extractOrthogonalPlanes.m*, *generateNIIFileName.m*, *readNII.m*, *visMotion.m*, *visReconstruction.m*, *visResiduals.m*, *visSegment.h*, *visTrajectory.m*, *writeData.m*, *writeNII.m*, *writeRaw.m*.

###### ./Libs
contains external MATLAB tools.

###### ./Lib/export_fig
from https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig

###### ./Lib/nextprod
from https://gist.github.com/fasiha/190203eac467d8b7f9ab2c83c3b3011e.

###### ./Lib/NIfTI_20140122
from https://uk.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

###### ./Lib/subtightplot
from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot

###### ./Methods
contains functions that implement generic methods for reconstruction: *aplGPU.m*, *blockGPU.m*, *build1DFTM.m*, *buildFilter.m*, *buildFoldM.m*, *buildHarmonicSampling.m*, *buildStandardDFTM.m*, *cdfFilt.m*, *compressCoils.m*, *compressMotion.m*, *computeROI.m*, *convertRotation.m*, *downsampleOperators.m*, *extractROI.m*, *dynInd.m*, *fctGPU.m*, *fct.m*, *fftGPU.m*, *filtering.m*, *flipping.m*, *fold.m*, *generateGrid.m*, *generateTransformGrids.m*, *groupwiseVolumeRegistration.m*, *ifctGPU.m*, *ifct.m*, *ifftGPU.m*, *ifold.m*, *ind2subV.m*, *indDim.m*, *mapMat.m*, *margosianFilter.m*, *matfun.m*, *mirroring.m*, *mirrorModulation.m*, *multDimMax.m*, *multDimMea.m*, *multDimMed.m*, *multDimMin.m*, *multDimSum.m*, *normm.m*, *numDims.m*, *parUnaFun.m*, *plugNoise.m*, *precomputeFactorsSincRigidTransform.m*, *pyramidPlan.m*, *removeOverencoding.m*, *resampling.m*, *resPop.m*, *resSub.m*, *restrictTransform.m*, *shearing.m*, *shifting.m*, *sincRigidTransformGradient.m*, *sincRigidTransform.m*, *sub2indV.m*, *upsampling.m*, *wrapToPiHalf.m*.

###### ./Reconstruction
contains functions to perform aligned reconstructions: *CGsolver.m*, *computeEnergy.m*, *constrain.m*, *CSsolver.m*, *decode.m*, *encodedecode.m*, *encode.m*, *LMsolver.m*, *precondition.m*, *regularize.m*, *solveC.m*, *solveG.m*, *solveXT.m*, *stopCondition.m*, *traceEstimation.m*.

###### ./Shearlets
contains functions that build the shearlet frames for regularization: *buildDirectionalFilter.m*, *buildQuadratureMirrorFilter.m*, *buildShearlet.m*, *buildWedgeBandpassAndLowpassFilters.m*, *checkFilterFeasibility.m*, *getShearletIdxs.m*.

###### ./Reconstruction
contains functions to simulate aligned reconstruction problems: *decodeDISORDER.m*, *directionDISORDER.m*, *encodeDISORDER.m*, *errorFit.m*, *generateParametersExp.m*, *gradientDISORDER.m*, *precondDISORDER.m*, *pyramidDISORDER.m*, *SNRToLevels.m*, *synthesizeEncoding.m*, *synthesizeT.m*, *synthesizeY.m*, *Xsolver.m*.


NOTE 1: Exemplary data is provided in the datasets *exampleSpectrum.txt* (used by *fig0304.m*), *xGT.mat* (used by *fig0506.m*), and *GT.mat* and *Q[1-3].mat* (used by *fig0708.m*). For runs without changing the paths, they should be placed in a folder
###### ../DISORDERData
Data generated when running the scripts appears in subfolders *Results5-6*, *Results7*. Nifti files are generated in subfolder *An-Ve* with suffixes *_Aq.nii* (no motion correction), *_Di.nii* (with correction) and *_Re.nii* (with robust correction).


NOTE 2: As for simulations, the provided execution should give a simplified version of the Figures in the paper, set 'quick=0' when calling *fig0506.m* to generate all the plots. As for reconstructions on real data *fig0708.m* these are computationally heavy; please, refer to the computation details in the manuscript to assess if your computing resources are adequate.


NOTE 3: The data to reproduce the experiments of Fig. 8 could not be linked to the release due to memory limitations. Please contact lucilio.cordero_grande@kcl.ac.uk if you are interested in these datasets (*BSSFP.mat*, *FLAIR.mat*, *MPRAGE.mat*, *SPGR.mat*, *TSE.mat*, used by *fig0708.m* with results in *Results8-A* and *Results8-B*).

