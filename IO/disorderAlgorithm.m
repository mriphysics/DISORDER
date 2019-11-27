function rec=disorderAlgorithm(rec)

%DISORDERALGORITHM   Sets the default parameters for DISORDER
%   REC=DISORDERALGORITHM(REC)
%   * REC is a reconstruction structure without the algorithmic parameters. 
%   At this stage it may contain the naming information (.Names), the 
%   status of the reconstruction (.Fail), the .lab information (.Par), the 
%   fixed plan information (.Plan) and the dynamic plan information (.Dyn) 
%   ** REC is a reconstruction structure with filled algorithmic 
%   information (.Alg)
%

%ALGORITHM (GENERIC)
if ~isfield(rec,'Alg');rec.Alg=[];end
if ~isfield(rec.Alg,'AlignedRec');rec.Alg.AlignedRec=2;end%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
if ~isfield(rec.Alg,'UseSoftMasking');rec.Alg.UseSoftMasking=1;end%To use soft instead of hard masking. If 2 we don't mask
if ~isfield(rec.Alg,'MargosianFilter');rec.Alg.MargosianFilter=1;end%To perform Margosian correction
if ~isfield(rec.Alg,'OverDec');rec.Alg.OverDec=ones(1,3);end%To overdecode the reconstructions
if ~isfield(rec.Alg,'WriteSnapshots');rec.Alg.WriteSnapshots=1;end%To write snapshots
if ~isfield(rec.Alg,'OnlyJSON');rec.Alg.OnlyJSON=0;end%To only write JSON metadata, not images
if ~isfield(rec.Alg,'EstimGFact');rec.Alg.EstimGFact=0;end%To estimate g-factors of reconstruction, not tested, use only experimentally

%DATA
if ~isfield(rec,'y');error('Data not present. Please provide data in field rec.y');end
if ~isfield(rec,'S');error('Sensitivities not present. Please provide sensitivities in field rec.S');end
if ~isfield(rec,'M');error('Mask not present. Please provide mask in the field rec.M');end

%OUTPUT FILE
if ~isfield(rec,'Names');rec.Names=[];end
if ~isfield(rec.Names,'pathOu') || isempty(rec.Names.pathOu);error('No output path provided. Please provide one in field rec.Names.pathOu');end
if ~isfield(rec.Names,'Name') || isempty(rec.Names.Name);error('No output file provided. Please provide one in field rec.Names.Name');end

%TRAJECTORIES
if ~isfield(rec,'Assign');rec.Assign=[];end
if ~isfield(rec.Assign,'z') || ~iscell(rec.Assign.z) || length(rec.Assign.z)<3 || isempty(rec.Assign.z{2}) || isempty(rec.Assign.z{3}) || numel(rec.Assign.z{2})~=numel(rec.Assign.z{3});error('Trajectories not set. These should go into a cell rec.Assign.z with second element for first PE (row vector) and third element for second PE (page vector). Both elements should have the same number of entries, which should correspond to the number of profiles');end

%SCAN INFO
if ~isfield(rec,'Par');rec.Par=[];end
if ~isfield(rec.Par,'Labels');rec.Par.Labels=[];end
if ~isfield(rec.Par.Labels,'TFEfactor')
    fprintf('Number of samples per effective shot (rec.Par.Labels.TFEfactor) not set. Default is used (only valid for steady-state)\n');
    rec.Par.Labels.TFEfactor=length(rec.Assign.z{2});
end    
if ~isfield(rec.Par.Labels,'ZReconLength')
    fprintf('Number of effective shots (rec.Par.Labels.ZReconLength) not set. Default is used (only valid for steady-state)\n');
    rec.Par.Labels.ZReconLength=1;
end

%RESOLUTION
if ~isfield(rec,'Enc');rec.Enc=[];end
if ~isfield(rec.Enc,'AcqVoxelSize')
    fprintf('Resolution (rec.Enc.AcqVoxelSize) not set. Default 1mm isotropic is used\n');
    rec.Enc.AcqVoxelSize=ones(1,3);%Resolution to 1mm
end
NY=size(rec.y);
if ~isfield(rec.Enc,'kRange');rec.Enc.kRange={[-NY(1) NY(1)-1],[-NY(2)/2 NY(2)/2-1],[-NY(3)/2 NY(3)/2-1]};end
if ~isfield(rec.Enc,'FOVSize');rec.Enc.FOVSize=NY(1:3);end
if ~isfield(rec.Enc,'AcqSize');rec.Enc.AcqSize=[2*NY(1) NY(2:3)];end

%GEOMETRY
if ~isfield(rec.Par,'Mine');rec.Par.Mine=[];end
if ~isfield(rec.Par.Mine,'Asca');rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);end
if ~isfield(rec.Par.Mine,'Arot');rec.Par.Mine.Arot=eye(4);end
if ~isfield(rec.Par.Mine,'Atra');rec.Par.Mine.Atra=eye(4);end
if ~isfield(rec.Par.Mine,'MTT');rec.Par.Mine.MTT=inv([-1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1]);end
if ~isfield(rec.Par.Mine,'APhiRec');rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;end

%SCAN GEOMETRY
if ~isfield(rec.Par,'Scan');rec.Par.Scan=[];end
if ~isfield(rec.Par.Scan,'MPS');rec.Par.Scan.MPS='HF LR PA';end
if ~isfield(rec.Par.Labels,'FatShiftDir');rec.Par.Labels.FatShiftDir='H';end
if ~isfield(rec.Par.Labels,'FoldOverDir');rec.Par.Labels.FoldOverDir='PA';end%This is the quick/first phase encode

%SCAN PARAMETERS
if ~isfield(rec.Par.Labels,'ScanDuration')
    fprintf('Scan duration not set. Default 1000s is used, temporal plots may not look well\n');
    rec.Par.Labels.ScanDuration=1000;
end
if ~isfield(rec.Par.Labels,'RepetitionTime')
    fprintf('Repetition time not set. Default 1s is used, temporal plots may not look well\n');    
    rec.Par.Labels.RepetitionTime=1;%1 second by default
end
if ~isfield(rec.Par.Labels,'TE')
    fprintf('Echo time not set. Default 1ms is used, temporal plots may not look well\n');    
    rec.Par.Labels.TE=0.001;%1 milisecond by default
end
if ~isfield(rec.Par.Labels,'FlipAngle');rec.Par.Labels.FlipAngle=0;end
if ~isfield(rec.Par.Scan,'Technique');rec.Par.Scan.Technique='T1TFE';end
if ~isfield(rec.Par.Mine,'Modal');rec.Par.Mine.Modal=7;end%Reconstruction modality

%ALGORITHM (SPECIFIC)
if ~isfield(rec.Alg,'parXT');rec.Alg.parXT=[];end
if ~isfield(rec.Alg.parXT,'NWend');rec.Alg.parXT.NWend=1;end%Number of within shot subdivisions at which to estimate motion, if multiple of 2 we assume tiling is used, otherwise rounded for temporal subdivisions
if ~isfield(rec.Alg.parXT,'accel');rec.Alg.parXT.accel=[1 0];end%Avoids motion estimation in the given number of levels of the spatio-temporal pyramid (to accelerate) - reconstruction in the given number of levels
if ~isfield(rec.Alg.parXT,'resolMax');rec.Alg.parXT.resolMax=4;end%2;%4;%Coarser resolution on which to compute motion (in mm)
if ~isfield(rec.Alg.parXT,'apod');rec.Alg.parXT.apod=[0.1 0 0];end%Apodization factors
if ~isfield(rec.Alg.parXT,'redFOV');rec.Alg.parXT.redFOV=1/3;end%Factor to reduce the FOV in the inferior direction for motion estimation
if ~isfield(rec.Alg.parXT,'perc');rec.Alg.parXT.perc=[0.9 0.9 0.9];end%Energy to be preserved after coil compression (if strictly lower than one, otherwise number of components preserved) / Redundant components used for motion estimation / Energy preserved when estimating motion
if ~isfield(rec.Alg.parXT,'traLimX');rec.Alg.parXT.traLimX=[0.05 0.02];end%Dispersion limit of translations/rotations for binning motion states
if ~isfield(rec.Alg.parXT,'traLimXT');rec.Alg.parXT.traLimXT=[0.05 0.02];end%Respectively translation and rotation limits for XT estimation
if ~isfield(rec.Alg.parXT,'echoesMin');rec.Alg.parXT.echoesMin=4;end%Minimum required number of samples to compute motion
if ~isfield(rec.Alg.parXT,'winit');rec.Alg.parXT.winit=1e-3;end%Initial weight of the quasi-Newton update
if ~isfield(rec.Alg.parXT,'meanT');rec.Alg.parXT.meanT=0;end%To constrain the parameters of the transforms solution to sum up to 0
if ~isfield(rec.Alg.parXT,'UseGiRingi');rec.Alg.parXT.UseGiRingi=1;end%To use Gibbs ringing filter
if ~isfield(rec.Alg.parXT,'GibbsRingi');rec.Alg.parXT.GibbsRingi=0.05;end%Gibbs ringing before SENSE unfolding
if ~isfield(rec.Alg.parXT,'discardHighRes');rec.Alg.parXT.discardHighRes=0.4;end%Those volumetric scans whose resolution is below this value (in mm) are not reconstructed
if ~isfield(rec.Alg.parXT,'convTransformJoint');rec.Alg.parXT.convTransformJoint=0;end%If 1 all transforms are estimated at each step, otherwise only those that have not converged
if ~isfield(rec.Alg.parXT,'tolerSolve');rec.Alg.parXT.tolerSolve=1e-3;end%Energy decrease in last iteration
if ~isfield(rec.Alg.parXT,'percRobustShot');rec.Alg.parXT.percRobustShot=[0.125 0.25];end%Percentiles for robust computation of expected inter-shot dispersion
if ~isfield(rec.Alg.parXT,'enerRobustShot');rec.Alg.parXT.enerRobustShot=0.95;end%Ratio of the error for acceptance of shots
if ~isfield(rec.Alg.parXT,'exploreMemory');rec.Alg.parXT.exploreMemory=0;end%To explore convergence without running the main methods
if ~isfield(rec.Alg.parXT,'writeInter');rec.Alg.parXT.writeInter=1;end%To write intermediate data
if ~isfield(rec.Alg.parXT,'computeCSRecon');rec.Alg.parXT.computeCSRecon=0;end%1;%To compute final CS-like reconstruction if different from zero. If not zero, it would set a resolution ratio of reconstructions with respect to baseline acquisition
if ~isfield(rec.Alg.parXT,'fractionOrder');rec.Alg.parXT.fractionOrder=0;end%Order for fractional derivative-based motion estimation
if ~isfield(rec.Alg.parXT,'saveFinal');rec.Alg.parXT.saveFinal=0;end%Save final results for further inspection
if ~isfield(rec.Alg.parXT,'groupSweeps');rec.Alg.parXT.groupSweeps=1;end%Factor to group sweeps for robust estimation
if ~isfield(rec.Alg.parXT,'maximumDynamics');rec.Alg.parXT.maximumDynamics=4;end%Maximum number of dynamics allowed
if ~isfield(rec.Alg,'parXB');rec.Alg.parXB=[];end
if ~isfield(rec.Alg.parXB,'dephaseCorrection');rec.Alg.parXB.dephaseCorrection=0;end%DCT radious for dephase correction, not in this version of the code
if ~isfield(rec.Alg.parXB,'tolerSolve');rec.Alg.parXB.tolerSolve=1;end%Minimum update of the phase per unit axial rotation (rad/rad), 1 radian in rotation producing 1 radian in phase would be acceptable
if ~isfield(rec.Alg.parXB,'winit');rec.Alg.parXB.winit=1e-3;end%Initial weight of the quasi-Newton update
if ~isfield(rec.Alg.parXB,'DCTdims');rec.Alg.parXB.DCTdims=[1 1 0];end%Initial weight of the quasi-Newton update


%DATA TYPE INFORMATION
if ~isfield(rec,'Plan');rec.Plan=[];end
if ~isfield(rec.Plan,'Dims');rec.Plan.Dims={'size','ky','kz','chan','dyn','card','echo','loca','mix','extr1','extr2','aver'};rec.Plan.NDims=length(rec.Plan.Dims);end
%1-> Raw / 2-> Empty / 3-> EPIRead / 4-> Navigator / 5-> Noise / 6-> Spectra
%7-> Sensitivities / 8-> Masks / 9-> G-facts / 10-> Chi2-maps / 11-> B0 / 12-> Reconstruction
%13-> Undistorted / 14-> Per-volume alignment / 15-> Per-excitation alignment / 16-> Volumetric alignment / 17-> Motion transform / 18-> Filtered reconstruction
%19-> Noise level / 20-> Number of components / 21-> Asymptotic error / 22-> Residuals / 23-> Frequency stabilization data / 24-> Image after frequency stabilization
%25-> Tracking information / 26->Tracked data / 27 -> SensitivityEigenMaps / 28 -> Chi2-maps in material coordinates
if ~isfield(rec.Plan,'Types');rec.Plan.Types={'z','','P','C','N','y', ...
            'S','M','G','X','B','x', ...
            'u','v','e','d','T','r', ...
            's','p','a','n','F','w', ...
            't','b','W','E'};end
if ~isfield(rec.Plan,'TypeNames');rec.Plan.TypeNames={'Ra','','Ny','Na','Ga','Sp', ...
                'Se','Ma','No','Ch','B0','Aq', ...
                'Un','Vo','Ex','Di','Tr','Re', ...
                'Si','Co','Ae','Er','Fr','Ws', ...
                'Fo','Vt','Ei','Ec'};rec.Plan.NTypes=length(rec.Plan.TypeNames);end
if ~isfield(rec.Plan,'CorrTypes');rec.Plan.CorrTypes={'random_phase','pda_fac','pda_index','meas_phase','sign'};rec.Plan.NCorrs=length(rec.Plan.CorrTypes);end

%ACTUALLY USED INFORMATION
if ~isfield(rec,'Dyn');rec.Dyn=[];end
if ~isfield(rec.Dyn,'Typ2Rec');rec.Dyn.Typ2Rec=[6;7;8];end%Types to reconstruct
if ~isfield(rec.Dyn,'Typ2Wri');rec.Dyn.Typ2Wri=zeros(1,28);rec.Dyn.Typ2Wri(12)=1;end%Types to write
if ~isfield(rec.Dyn,'Debug');rec.Dyn.Debug=2;end%To show debug information

%ERROR HANDLING
if ~isfield(rec,'Fail');rec.Fail=0;end%Flag to abort

%GPU
if ~isfield(rec.Dyn,'BlockGPU');rec.Dyn.BlockGPU=1-double(gpuDeviceCount>0);end%Flag to block GPU computations
if ~isfield(rec.Dyn,'GPU');rec.Dyn.GPU=2;end%Usage of GPU, 1 if arrays already in GPU, 2 to convert (less more memory demanding), 0 for not using it
if ~isfield(rec.Dyn,'MaxMem');rec.Dyn.MaxMem=[6e6 2e6 1e6];end%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing

%OTHER
if ~isfield(rec.Plan,'SuffOu');rec.Plan.SuffOu='';end%Suffix to add to output data
if ~isfield(rec.Plan,'Suff');rec.Plan.Suff='';end%Suffix to add to images
if ~isfield(rec.Par.Mine,'Proce');rec.Par.Mine.Proce=0;end%Flag to indicate that this is a reconstruction problem instead of a preprocessing problem
if ~isfield(rec.Par.Mine,'Signs');rec.Par.Mine.Signs=[];end%For rotating PEs, not used
if ~isfield(rec.Par.Mine,'Nat');rec.Par.Mine.Nat=1;end%Native PE in case of rotating PEs, not used
if ~isfield(rec.Par,'Encoding');rec.Par.Encoding=[];end%Additional encoding information, not used

