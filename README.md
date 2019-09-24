# DISORDER
Tools for aligned reconstructions of volumetric MRI data

This repository provides tools to implement the methods and reproduce the experiments included in the manuscript ''Motion corrected MRI with DISORDER: Distributed and Incoherent Sample Orders for Reconstruction Deblurring using Encoding Redundancy'', L Cordero-Grande, G Ferrazzi, RPAG Teixeira, J O'Muircheartaigh, AN Price, and JV Hajnal, arXiv, 2019.

The code has been developed in MATLAB and has the following structure:

###### ./
contains the functions to run the illustrations and experiments included in Figs. 3-9  of the manuscript: *fig[0304,050607,0809].m*.

###### ./Acquisition
contains template functions to implement the DISORDER reorderings in MRI scanners: *electrostaticRepulsionDISORDER.m*, *samplingDISORDER.m*.

###### ./Libs
contains external MATLAB tools.

###### ./Lib/nextprod
from https://gist.github.com/fasiha/190203eac467d8b7f9ab2c83c3b3011e.

###### ./Lib/NIfTI_20140122
from https://uk.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

###### ./Lib/subtightplot
from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot

###### ./Methods
contains functions that implement generic methods for reconstruction: *aplGPU.m*, *blockGPU.m*, *build1DFTM.m*, *buildFilter.m*, *buildFoldM.m*, *buildHarmonicSampling.m*, *buildStandardDFTM.m*, *cdfFilt.m*, *compressCoils.m*, *compressMotion.m*, *convertRotation.m*, *downsampleOperators.m*, *dynInd.m*, *fctGPU.m*, *fct.m*, *fftGPU.m*, *filtering.m*, *flipping.m*, *fold.m*, *generateGrid.m*, *generateTransformGrids.m*, *groupwiseVolumeRegistration.m*, *ifctGPU.m*, *ifct.m*, *ifftGPU.m*, *ifold.m*, *ind2subV.m*, *indDim.m*, *mapMat.m*, *margosianFilter.m*, *matfun.m*, *mirroring.m*, *mirrorModulation.m*, *multDimMax.m*, *multDimMea.m*, *multDimMed.m*, *multDimSum.m*, *normm.m*, *numDims.m*, *parUnaFun.m*, *plugNoise.m*, *precomputeFactorsSincRigidTransform.m*, *pyramidPlan.m*, *removeOverencoding.m*, *resampling.m*, *resPop.m*, *resSub.m*, *restrictTransform.m*, *shearing.m*, *shifting.m*, *sincRigidTransformGradient.m*, *sincRigidTransform.m*, *sub2indV.m*, *upsampling.m*, *wrapToPiHalf.m*.

###### ./Reconstruction
contains functions to perform aligned reconstructions: *CGsolver.m*, *computeEnergy.m*, *constrain.m*, *CSsolver.m*, *decode.m*, *encodedecode.m*, *encode.m*, *LMsolver.m*, *precondition.m*, *readNII.m*, *regularize.m*, *solveC.m*, *solveG.m*, *solveXT.m*, *stopCondition.m*, *traceEstimation.m*, *writeData.m*, *writeNII.m*.

###### ./Shearlets
contains functions that build the shearlet frames for regularization: *buildDirectionalFilter.m*, *buildQuadratureMirrorFilter.m*, *buildShearlet.m*, *buildWedgeBandpassAndLowpassFilters.m*, *checkFilterFeasibility.m*, *getShearletIdxs.m*.

###### ./Reconstruction
contains functions to simulate aligned reconstruction problems: *decodeDISORDER.m*, *directionDISORDER.m*, *encodeDISORDER.m*, *errorFit.m*, *generateParametersExp.m*, *gradientDISORDER.m*, *precondDISORDER.m*, *pyramidDISORDER.m*, *SNRToLevels.m*, *synthesizeEncoding.m*, *synthesizeT.m*, *synthesizeY.m*, *Xsolver.m*.


NOTE 1: Exemplary data is provided in the datasets *BSSFP.mat*,*exampleSpectrum.txt*, *FLAIR.mat*, *GT.mat*, *MPRAGE.mat*, *Q[1-3].mat*, *SPGR.mat*, *TSE.mat*, *xGT.mat*. For runs without changing the paths, they should be placed in a folder
###### ../DISORDERData
Data generated when running the scripts is also stored in this folder as *BSSFP_[Aq,AqPh,Di,DiPh,Re,RePh].nii*, *fig0[5-7].mat*, *FLAIR_[Aq,AqPh,Di,DiPh,Re,RePh].nii*, *GT_[Aq,AqPh,Di,DiPh,Re,RePh].nii*, *MPRAGE_[Aq,AqPh,Di,DiPh,Re,RePh].nii*, *Q[1-3]_[Aq,AqPh,Di,DiPh,Re,RePh].nii*, *SPGR_[Aq,AqPh,Di,DiPh,Re,RePh].nii*, *TSE_[Aq,AqPh,Di,DiPh,Re,RePh].nii*.


NOTE 2: As for simulations, the provided execution should give a simplified version of the Figures in the paper, set 'quick=0' when calling *fig050607.m* to generate all the plots. As for reconstructions on real data *fig0809.m* these are computationally heavy; please, refer to the computation details in the manuscript to assess if your computing resources are adequate.

