function rec=solveXT(rec)

%SOLVEXT   Performs an aligned reconstruction according to the manuscript
%"Motion corrected magnetic resonance imaging with DISORDER: Distributed 
%and Incoherent Sample Orders for Reconstruction Deblurring using Encoding 
%Redundancy", L Cordero-Grande, G Ferrazzi, RPAG Teixeira, AN Price, 
%J O’Muircheartaigh, and JV Hajnal, 2019.
%   REC=SOLVEXT(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   ** REC is a reconstruction structure with reconstructed data and
%   surrogate information (rec.(rec.Plan.Types))
%

tsta=tic;

%SET DEFAULT PARAMETERS
rec=disorderAlgorithm(rec);

%FILENAME
folderSnapshots=strcat(rec.Names.pathOu,filesep,'An-Ve_Sn');
[~,fileName]=fileparts(generateNIIFileName(rec));

%SHORTCUTS FOR COMMONLY USED STRUCTURE FIELDS AND INITIALIZERS
voxSiz=rec.Enc.AcqVoxelSize;%Acquired voxel size
parXT=rec.Alg.parXT;
gpu=rec.Dyn.GPU;gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);if gpuIn;gpuFIn=2;else gpuFIn=0;end
if size(rec.y,5)>=parXT.maximumDynamics;fprintf('Problem too big to fit in memory (%d dynamics for a limit of %d)\n',size(rec.y,5),parXT.maximumDynamics);rec.Fail=1;return;end
on=cell(1,4);for n=1:4;on{n}=ones(1,n);end
typ2Rec=rec.Dyn.Typ2Rec;

%RECONSTRUCTION PLAN
[resPyr,L,estT,resIso]=pyramidPlan(voxSiz,parXT.resolMax,parXT.NWend,parXT.accel);

%ROI COMPUTATION AND EXTRACTION IN THE READOUT DIRECTION
rec.Enc.ROI=computeROI(rec.M,[],[0 0 0],[1 0 0]);
if rec.Dyn.Debug>=2
    fprintf('ROI:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
    fprintf('Acquisition directions: %s\n',rec.Par.Scan.MPS);
end
for n=typ2Rec';datTyp=rec.Plan.Types{n};
    if ~ismember(n,[5 12]);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,1);end
end

%PERMUTE PHASE ENCODES TO THE FIRST DIMENSIONS
perm=1:rec.Plan.NDims;perm(1:3)=[3 2 1];
for n=typ2Rec'        
    if ~ismember(n,[5 12]);datTyp=rec.Plan.Types{n};
        rec.(datTyp)=permute(rec.(datTyp),perm);
    end
end
voxSiz=voxSiz(perm(1:3));parXT.apod=parXT.apod(perm(1:3));

%COIL ARRAY COMPRESSION AND RECONSTRUCTED DATA GENERATION
NX=size(rec.M);
[S,y,eivaS]=compressCoils(rec.S,parXT.perc,rec.y);
NS=size(S);
if gpuIn;y=gpuArray(y);end
if rec.Dyn.Debug>=2 && ~isempty(parXT.perc)
    if parXT.perc(1)<1;fprintf('Number of compressed coil elements at%s%%: %d (%s )\n',sprintf(' %0.2f',parXT.perc*100),NS(4),sprintf(' %d',eivaS));else fprintf('Number of compressed coil elements: %d (%s )\n',NS(4),sprintf(' %d',eivaS));end
end

%APODIZE + REARRANGE DATA + CORRECT FOR INTERLEAVED REPEATS
NY=size(y);NY(end+1:rec.Plan.NDims)=1;
y=bsxfun(@times,y,ifftshift(buildFilter(NY(1:3),'tukey',[10 10 1],gpuIn,parXT.apod)));%Apodize
[y,NY]=resSub(y,5:rec.Plan.NDims);NY(end+1:5)=1;%Different averages along the 5th dimension
y=gather(y);

%ACQUISITION STRUCTURE
NEchos=max(rec.Par.Labels.TFEfactor,rec.Par.Labels.ZReconLength);
NProfs=numel(rec.Assign.z{2});
if NProfs<=1 || NProfs~=numel(rec.Assign.z{3});fprintf('SEPARABLE Y-Z TRAJECTORIES. ALIGNED RECONSTRUCTION IS NOT PERFORMED\n');rec.Fail=1;return;end
if sum(cdfFilt(abs(diff(rec.Assign.z{2}(:))),'med'))<sum(cdfFilt(abs(diff(rec.Assign.z{3}(:))),'med'));fprintf('Slow direction is: %s\n',rec.Par.Scan.MPS(4:5));else fprintf('Slow direction is: %s\n',rec.Par.Scan.MPS(7:8));end
kTraj=zeros([NProfs 2],'single');
for n=1:2;kTraj(:,n)=rec.Assign.z{perm(n)}(:);end
if NEchos==1;NEchos=NProfs;end
NShots=NProfs/NEchos;
if mod(NShots,1)~=0;fprintf('NUMBER OF SHOTS: %.2f, PROBABLY AN INTERRUPTED SCAN\n',NShots);rec.Fail=1;return;end
NRepeats=NY(5);
if NShots<NRepeats
    NShots=NRepeats;
    NProfs=NShots*NEchos;
    kTraj=repmat(kTraj,[NRepeats 1]);
    NEchos=NProfs/NShots;
end
if rec.Dyn.Debug>=2
    fprintf('Number of echoes: %d\n',NEchos);
    fprintf('Number of shots: %d\n',NShots);
end
kTrajSS=reshape(kTraj,[NEchos NShots 2]);
%PLOT TRAJECTORIES
if rec.Alg.WriteSnapshots;visTrajectory(kTrajSS,2,strcat(folderSnapshots,filesep,'Trajectories'),strcat(fileName,rec.Plan.SuffOu));end
kRange=single(zeros(2,2));
for n=1:2;kRange(n,:)=rec.Enc.kRange{perm(n)};end
kShift=(floor((diff(kRange,1,2)+1)/2)+1)';
kSize=(diff(kRange,1,2)+1)';
kIndex=bsxfun(@plus,kTraj,kShift);
rec.Par.Mine.DISORDER.kTraj=gather(kTraj);
rec.Par.Mine.DISORDER.kRange=gather(kRange);
rec.Par.Mine.DISORDER.kShift=gather(kShift);
rec.Par.Mine.DISORDER.kSize=gather(kSize);
rec.Par.Mine.DISORDER.kIndex=gather(kIndex);
%STEADY-STATE CORRECTIONS
isSteadyState=0;
if mod(NShots,NRepeats)==0
    NShotsUnique=NShots/NRepeats;  
    isSteadyState=(NShotsUnique==1);%No shots have been identified
end
fprintf('Steady-state: %d\n',isSteadyState);
if isSteadyState
    kTrajs=abs(fct(kTraj));
    N=size(kTrajs);
    [~,iAs]=max(dynInd(kTrajs,1:floor(N(1)/2),1),[],1);
    iM=min(iAs,[],2);
    NSamples=NProfs/(iM-1);
    NSweeps=round(NProfs/(NSamples*parXT.groupSweeps));
    NSamples=NProfs/NSweeps;
    fprintf('Number of sweeps: %d\n',NSweeps);
else
    NSamples=NEchos;
    NSweeps=NShots;
end
NStates=NSweeps;
%SWEEP SUBDIVISION AND TIME INDEX
sweepSample=ceil(NSweeps*(((1:NProfs)-0.5)/NProfs));
stateSample=sweepSample;
timeIndex=zeros([NY(1:2) NY(5)]);
if gpu;timeIndex=gpuArray(timeIndex);end
for n=1:size(kIndex,1)
    for s=1:NY(5)
        if timeIndex(kIndex(n,1),kIndex(n,2),s)==0
            timeIndex(kIndex(n,1),kIndex(n,2),s)=n;
            break;
        end
    end
end
for m=1:2;timeIndex=ifftshift(timeIndex,m);end
if isSteadyState
    timeSample=(0:NProfs-1)*rec.Par.Labels.RepetitionTime(1)/1000;
else
    TPerShot=rec.Par.Labels.ScanDuration/NShots;    
    timeSample=(sweepSample-1)*TPerShot;
    if strcmp(rec.Par.Scan.Technique,'TSE') || strcmp(rec.Par.Scan.Technique,'TIR')
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*(rec.Par.Labels.TE(1)/(1000*NEchos/2));
    else
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*rec.Par.Labels.RepetitionTime(1)/1000;
    end
end
sweepToSample=cell(1,NSweeps);
sampleToSweep=1:NSweeps;
for n=1:NSweeps;sweepToSample{n}=n;end

%CREATE RECONSTRUCTION DATA ARRAY
NX=size(rec.M);NX(end+1:3)=1;
rec.x=zeros(NX,'single');%Uncorrected non-regularized
rec.d=zeros([NX L-sum(resPyr~=1)],'single');%Corrected non-regularized
if parXT.computeCSRecon;rec.r=zeros([NX L-sum(resPyr~=1)],'single');end%Corrected regularized
x=zeros(NX,'single');if gpuIn;x=gpuArray(x);end
typ2RecI=[12;16];
if parXT.computeCSRecon;typ2RecI=[typ2RecI;18];end
for n=typ2RecI'
    if ~any(rec.Dyn.Typ2Rec==n);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,n);rec.Dyn.Typ2Wri(n)=1;end
end
typ2Rec=rec.Dyn.Typ2Rec;

%INITIALIZE TRANSFORM, NUMBER OF EXTERN ITERATIONS, TOLERANCES, IDENTIFY 
%STEADY STATE SEQUENCES 
rec.T=zeros([on{4} NSweeps 6],'single');
subT=ones(1,NStates);%To perform subdivisions of T
nExtern=99;
tolType={'Energy','RelativeResidual2Error'};%For CG--without vs with motion correction
tol=[parXT.tolerSolve 0];%For CG
nIt=[300 1];%For CG
nowithin=0;%To stop estimation of within motion

%ENERGY
En=cell(1,L);Re=cell(1,L);EnX=[];
mSt=cell(1,L);iSt=cell(1,L);
tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time arranging: %.3f s\n',tend);end

%SOLVE WITHOUT MOTION
tsta=tic;
resAni=single(on{3});
[~,indMinRes]=min(voxSiz);        
resAni(indMinRes)=voxSiz(indMinRes)/resPyr(end);
indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes));   
resAni=voxSiz./resAni;
BlSz=max(sqrt(prod(voxSiz(1:2)./resAni(1:2))),1);%Block sizes baseline
E.bS=round([2 4].^log2(BlSz));E.oS=on{2};%Block sizes for our GPU
if rec.Dyn.Debug>=2;fprintf('Block sizes:%s\n',sprintf(' %d',E.bS));end
E.Je=1;%To use joint encoding/decoding
E.Uf=cell(1,3);for n=1:3;E.Uf{n}.NY=NY(n);E.Uf{n}.NX=NX(n);end%Folding
EH.Ub=cell(1,3);for n=1:3;EH.Ub{n}.NY=NY(n);EH.Ub{n}.NX=NX(n);end%Unfolding
[E.UAf,EH.UAb]=buildFoldM(NX(1:3),NY(1:3),gpuIn,1);%Folding/Unfolding matrix (empty is quicker!)
R=[];%Regularizer
ncx=1:eivaS(1);

%MASKING
CX.Ma=rec.M;
if gpuIn;CX.Ma=gpuArray(CX.Ma);end
if rec.Alg.UseSoftMasking==2;CX.Ma(:)=1;end
discMa=all(ismember(CX.Ma(:),[0 1]));
if ~discMa%Soft masking
    fprintf('Using soft masking\n');
    M=buildFilter(2*[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX],'tukeyIso',1,gpuIn,1,1);
    CX.Ma=abs(filtering(CX.Ma,M,1));M=[];
    if size(S,6)>1 || isfield(rec,'W')%We do not unfold in the PE direction
        NMa=[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX];
        for m=1:3
            if ~isempty(E.Uf{m});CX.Ma=fold(CX.Ma,n,E.Uf{m}.NX,E.Uf{m}.NY,E.UAf{m});end
        end         
        CX.Ma=resampling(CX.Ma,NMa,2);
    end
    if isfield(rec,'W');CX.Ma=bsxfun(@times,CX.Ma,rec.W);end
    Ti.la=CX.Ma;
    R.Ti.la=1./(abs(Ti.la).^2+0.01);    
    CX=[];
end

%SENSITIVITIES
E.Sf=dynInd(S,ncx,4);E.dS=[1 ncx(end)];
if gpuIn;E.Sf=gpuArray(E.Sf);end

%PRECONDITION
if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
else P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
end

yX=mean(dynInd(y,ncx,4),5);
if gpuIn;yX=gpuArray(yX);end
gibbsRing=parXT.UseGiRingi*parXT.GibbsRingi;
if gibbsRing~=0 && gibbsRing<=1;yX=filtering(yX,buildFilter(NY(1:3),'tukeyIso',[],gpuIn,gibbsRing));end
nX=nIt(1);
if parXT.exploreMemory;nX=0;end%To test flow without running the main methods       

if nX~=0
    [xRes,EnX]=CGsolver(yX,E,EH,P,CX,R,x,nX,tolType{1},tol(1));
    rec.x=gather(xRes);
end
if ~isempty(EnX)%Print and store final energy
    EnX=EnX/numel(yX);
    if rec.Dyn.Debug>=2;fprintf('Energy evolution last step:%s\n',sprintf(' %0.6g',sum(EnX,1)));end  
    EnX=[];
end;yX=[];

if gpuIn;y=gpuArray(y);end
Residuals=encode(rec.x,E)-y;
for m=1:2;Residuals=fftGPU(Residuals,m)/sqrt(NY(m));end
Residuals=normm(Residuals,[],3:4);
if rec.Alg.WriteSnapshots
    if strcmp(rec.Par.Labels.FatShiftDir,'F');E.nF=1:floor((1-parXT.redFOV)*NX(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.nF=floor(1+parXT.redFOV*NX(3)):NX(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NX(3);end
    visSegment(flip(permute(dynInd(rec.x,E.nF,3),[2 1 3]),2),[],2,1,[],[],'Uncorrected',strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.SuffOu,'_Aq'));
    visResiduals([],[],[],Residuals,kTraj,2,{strcat(folderSnapshots,filesep,'ResidualsStates'),strcat(folderSnapshots,filesep,'ResidualsSpectrum')},strcat(fileName,rec.Plan.SuffOu,sprintf('_l=%d',0)));
end
E.Sf=[];

tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing non-corrected reconstruction: %.3f s\n',tend);end

%ITERATING THROUGH THE DIFFERENT LEVELS
for l=1:L    
    %WE RESET MOTION CONVERGENCE INSIDE
    nReset=0;THist=[];
    leff=l-sum(resPyr~=1);

    %UPSAMPLE THE TRANSFORM AND GENERATE MOTION PARAMETERS   
    rec.T=repelem(rec.T,1,1,1,1,subT,1);
    subTT=repelem(subT,1,subT);
    sampleToSweep=repelem(sampleToSweep,1,subT);
    sweepToSample=cell(1,NSweeps);
    c=0;
    for n=1:length(subT)
        if subT(n)==2
            indState=find(stateSample==(n+c));
            Ni=length(indState);
            indState=indState(floor(Ni/2)+1);
            stateSample(indState:end)=stateSample(indState:end)+1;
            sweepToSample{sampleToSweep(n+c)}=[sweepToSample{sampleToSweep(n+c)} n+c:n+c+1];
            c=c+1;
        else
            sweepToSample{sampleToSweep(n+c)}=[sweepToSample{sampleToSweep(n+c)} n+c];
        end
    end 
    timeState=regionprops(stateSample,timeSample,'MeanIntensity');
    timeStateLimitInf=regionprops(stateSample,timeSample,'MinIntensity');
    timeStateLimitSup=regionprops(stateSample,timeSample,'MaxIntensity');
    timeState={timeState.MeanIntensity};
    timeStateLimitInf={timeStateLimitInf.MinIntensity};
    timeStateLimitSup={timeStateLimitSup.MaxIntensity};
    timeState=cat(1,timeState{:});
    timeStateLimitInf=cat(1,timeStateLimitInf{:});
    timeStateLimitSup=cat(1,timeStateLimitSup{:});
    if leff<=1;timeSweep=timeState;end
    %NOTE THIS ONLY WORKS IF THERE IS NO INTERSHOT MOTION
    if leff<=1
        if ~isSteadyState;timeState=cat(3,timeStateLimitInf,timeStateLimitSup);else timeState=cat(3,timeStateLimitInf,timeState,timeStateLimitSup);end
    end
    
    %RELATIVE SPATIAL RESOLUTION CALCULATIONS AT THIS LEVEL
    resAni=single(on{3});
    [~,indMinRes]=min(voxSiz);        
    resAni(indMinRes)=voxSiz(indMinRes)/resPyr(l);
    indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
    resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes));   
    resAni=voxSiz./resAni;    

    NT=size(rec.T);    
    NStates=NT(5);
    subT=ones(1,NStates);%To perform subdivisions of T
    cT=zeros(NT(1:5),'single');%Flag to indicate whether a given motion state can be assumed to have converged     
    w=parXT.winit*ones(NT(1:5),'single');%Weights of LM iteration            
    NSamplesPerState=accumarray(stateSample(:),ones(NProfs,1))';    
    if rec.Dyn.Debug>=2
        fprintf('-----------------------\n');
        fprintf('Resolution level: %d\n',l);
        fprintf('Number of motion states: %d\n',NStates);
        fprintf('Minimum/Maximum number of echoes: %d/%d\n',min(NSamplesPerState),max(NSamplesPerState));
        if any(NSamplesPerState(:)<parXT.echoesMin);fprintf('Trying to operate below the minimum required number of reads, motion estimates will not be obtained for stability\n');end
    end
    if estT(l)
        cT(NSamplesPerState(:)<parXT.echoesMin)=1;
        cT=reshape(cT,NT(1:5));
    end

    %ENCODING STRUCTURES AND BLOCK SIZES FOR GPU COMPUTATION AT THIS LEVEL
    if prod(round(NX(1:3).*resAni))<=rec.Dyn.MaxMem(2) && gpuIn;S=gpuArray(S);rec.M=gpuArray(rec.M);end
    if ~discMa
        [xRes,SRes,yRes,timeIndexRes,R.Ti.la,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,S,y,timeIndex,Ti.la,gibbsRing,0,gpuIn);
        R.Ti.la=1./(abs(R.Ti.la).^2+0.01);
    else
        [xRes,SRes,yRes,timeIndexRes,CX.Ma,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,S,y,timeIndex,rec.M,gibbsRing,0,gpuIn);
    end

    NXres=size(xRes);NXres(end+1:3)=1;NYres=size(yRes);NYres(end+1:5)=1;
    if prod(round(NX(1:3).*resAni))<=rec.Dyn.MaxMem(2) && gpuIn;S=gather(S);rec.M=gather(rec.M);end
    if l==L;S=[];y=[];end
    for n=1:2;yRes=fftGPU(yRes,n)/sqrt(NYres(n));end

    %WEIGHTING FUNCTION FOR MOTION ESTIMATION (ROI)
    if strcmp(rec.Par.Labels.FatShiftDir,'F');E.nF=1:floor((1-parXT.redFOV)*NXres(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.nF=floor(1+parXT.redFOV*NXres(3)):NXres(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NXres(3);end
    EH.Mbe=zeros([1 1 NYres(3)],'like',real(xRes));EH.Mbe(E.nF)=1;%To use the extracted region for energy computation

    %HARMONIC SAMPLING
    E.mSt=timeIndexRes;  
    E.mSt(timeIndexRes~=0)=stateSample(timeIndexRes(timeIndexRes~=0));

    if estT(l) || ~isfield(EH,'We');EH.We=[];end    
    if parXT.fractionOrder~=0;GRes=buildFilter([NXres(1:2) length(E.nF)],'FractionalFiniteDiscreteIso',NX./NXres,gpuIn,parXT.fractionOrder);else GRes=[];end%For fractional finite difference motion estimation           
    [ySt,FSt,GRes,E.iSt,E.nSt,E.mSt,nSa]=buildHarmonicSampling(yRes,GRes,E.mSt,NStates);iSt{l}=cat(1,E.iSt{:});yRes=gather(yRes);%Samples as cell arrays
    if any(nSa(1:NStates)==0);nowithin=1;fprintf('Probably not DISORDER, returning without performing corrections\n');return;end       

    if leff<=1;mStSweeps=E.mSt;end
    BlSz=max(sqrt(prod(voxSiz(1:2)./resAni(1:2))),1);%Block sizes baseline
    E.bS=round([2 4].^log2(BlSz));E.oS=on{2};%Block sizes for our GPU
    if rec.Dyn.Debug>=2;fprintf('Block sizes:%s\n',sprintf(' %d',E.bS));end
    E.Je=1;%To use joint encoding/decoding
    E.Uf=cell(1,3);for n=1:3;E.Uf{n}.NY=NYres(n);E.Uf{n}.NX=NXres(n);end%Folding
    EH.Ub=cell(1,3);for n=1:3;EH.Ub{n}.NY=NYres(n);EH.Ub{n}.NX=NXres(n);end%Unfolding
    [E.UAf,EH.UAb]=buildFoldM(NXres(1:3),NYres(1:3),gpuIn,1);%Folding/Unfolding matrix (empty is quicker!)

    %CONTROL OF TRANSFORM CONVERGENCE / CONSTRAIN THE TRANSFORM / ROI IN
    %READOUT FOR MOTION ESTIMATION
    E.w=w;E.cT=cT;E.winit=parXT.winit;
    if parXT.meanT;CT.mV=5;else CT=[];end
    E.tL=parXT.traLimXT*resIso(l);%Minimum update to keep on doing motion estimation

    cont=0; 
    if isfield(E,'En');E=rmfield(E,'En');end    

    %ITERATIVE SOLVER
    for n=1:nExtern
        if rec.Dyn.Debug>=2;fprintf('Iteration: %d\n',n);end
        %RESET CONVERGENCE OF MOTION
        if mod(cont,nReset)==0 || all(E.cT(:)) || parXT.convTransformJoint
            nReset=nReset+1;cont=0;%We increment the number of iterations till next reset
            E.cT(:)=0;%We reset
            E.cT(cT==1)=1;%We fix to 1 the low frequency states in case that for whatever reason cT was set to 1
            if rec.Dyn.Debug>=2;fprintf('Resetting motion convergence\n');end
        end
        if (n==1 || l<L) && any(subTT(:)==2)           
            E.cT(subTT==1)=1;
        end%For within-shot in the first iteration we only estimate those coming from outliers
        if rec.Dyn.Debug>=2;fprintf('Explored motion states: %d of %d\n',NStates-sum(single(E.cT)),NStates);end
        cont=cont+1;

        %TRANSFORM PARAMETERS PER SEGMENT/OUTLIERS/DATA/HARMONICS/SEGMENT INDEXES
        [E.Tr,E.Fs,E.nEc]=compressMotion(rec.T,FSt,resIso(end),parXT);
        if rec.Dyn.Debug>=2;fprintf('Number of binned motion states: %d of %d\n',size(E.Tr,5),NT(5));end

        %SEGMENTS / MOTION STATES / BLOCK SIZES / TRANSFORM FACTORS 
        if ~isempty(E.Fs);E.NSe=length(E.Fs{1});E.NMs=size(E.Tr,5);else E.NSe=1;E.NMs=1;end%Number of segments (including out of elliptic shutter) / Number of motion states            
        E.NRe=E.NSe-E.NMs;%Number of outern segments
        E.dS(1)=E.NSe;

        %CG SOLVER
        ncx=1:eivaS(1+2*estT(l));
        E.Sf=dynInd(SRes,ncx,4);E.dS(2)=ncx(end);
        if gpuIn;E.Sf=gpuArray(E.Sf);end

        %PRECONDITION
        if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
        else P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
        end            

        yX=dynInd(ySt,ncx,4);            
        %LAST ITERATION OUTLIER REJECTION
        if (leff>0 && estT(l)==0) && rec.Alg.AlignedRec==2 && exist('We','var') 
            EH.We=ones([size(yX,1) 1],'like',Residuals);
            for s=1:NSweeps;EH.We(mStSweeps==s)=We(s);end%Soft weights                
        end

        if estT(l)
            if l==1 && n==1;nX=nIt(2);
            elseif l==1;nX=max(nX-1,1);
            elseif resPyr(l)~=resPyr(l-1) && n==1;nX=max(nIt(2)-1,1);
            elseif resPyr(l)~=resPyr(l-1);nX=max(nX-1,1);
            elseif resPyr(l)==resPyr(l-1) && n==1;nX=0;
            elseif resPyr(l)==resPyr(l-1) && n==2;nX=max(nIt(2)-2,1);
            elseif resPyr(l)==resPyr(l-1);nX=max(nX-(n-2),1);
            end
        else
            nX=nIt(1);
        end
        if parXT.exploreMemory;nX=0;end%To test flow without running the main methods       
        %if nX~=0 && (leff<=0 || n==1);[xRes,EnX]=CGsolver(yX,E,EH,P,CX,R,xRes,nX,tolType{estT(l)+1},tol(estT(l)+1));end
        if nX~=0;[xRes,EnX]=CGsolver(yX,E,EH,P,CX,R,xRes,nX,tolType{estT(l)+1},tol(estT(l)+1));end

        if ~isempty(EnX)%Print and store final energy
            EnX=EnX/numel(yX);
            if rec.Dyn.Debug>=2;fprintf('Energy evolution last step:%s\n',sprintf(' %0.6g',sum(EnX,1)));end  
            EnX=[];
        end

        %ENERGY COMPUTATION
        E.Sf=dynInd(SRes,ncx,4);E.dS(2)=ncx(end);
        if gpuIn;E.Sf=gpuArray(E.Sf);end
        yX=dynInd(ySt,ncx,4);
        [Re{l}(:,n),~,fullResid]=computeEnergy(yX,xRes,E,[],EH,[],[],[],[],1);Re{l}(:,n)=Re{l}(:,n)/numel(yX);
        En{l}(:,n)=normm(yX,[],setdiff(1:numDims(yX),1))/numel(yX);
        if rec.Dyn.Debug>=2;fprintf('Residuals: %0.6g\n',sum(Re{l}(:,n)));end        

        %ENERGY PER STATES
        fullResid=gather(dynInd(fullResid,iSt{l},1,fullResid));
        Residuals=Re{l}(:,n);Residuals(iSt{l})=Residuals;
        Residuals=reshape(Residuals,NYres([1:2 5]));
        mStSort=mStSweeps;mStSort(iSt{l})=mStSort;
        Sweeps=reshape(mStSort,NYres([1:2 5]));         

        %OUTLIER DETECTION
        if leff>0;percUse=(0.15:0.05:0.85)*100;
        else percUse=(0.3:0.05:0.7)*100;
        end
        WeO=zeros([NSweeps length(percUse)],'like',Residuals);
        for s=1:NSweeps;WeO(s,:)=prctile(log(Residuals(Sweeps==s)),percUse);end
        WeO=permute(WeO,[1 3 2]);           
        Westd=diff(prctile(WeO,parXT.percRobustShot*100,1),1,1)./diff(norminv(parXT.percRobustShot)); 
        Wmea=prctile(WeO,mean(parXT.percRobustShot)*100,1);
        Wmea=Wmea+Westd*sqrt(2)*erfcinv(2*mean(parXT.percRobustShot));%sqrt(2)*erfcinv(2*mean(parXT.percRobustShot))=norminv(mean(parXT.percRobustShot))
        We=(WeO-Wmea);
        We=We./Westd;
        We=mean(We,3);         
        We=min((1-normcdf(We))/((1-parXT.enerRobustShot)/NSweeps),1);     
        outlWe=We<1;

        if rec.Alg.WriteSnapshots && estT(l)==0;visResiduals(We,outlWe,timeSweep,resampling(Residuals,NY(1:2),1),kTraj,2,{strcat(folderSnapshots,filesep,'ResidualsStates'),strcat(folderSnapshots,filesep,'ResidualsSpectrum')},strcat(fileName,rec.Plan.SuffOu,sprintf('_l=%d',l)));end
        %CS-LIKE/ROBUST RECONSTRUCTION
        if leff>0 && estT(l)==0
            subT(:)=1;
            if rec.Alg.AlignedRec<4;outlWeV=find(outlWe);else outlWeV=1:length(outlWe);end%AlignedRec==4 explore all, AlignedRec==3, explore only artifacted
            for o=1:length(outlWeV);subT(sweepToSample{outlWeV(o)})=2;end

            if all(subT(:)==1) || nowithin;L=l;end
        end
        if (leff>0 && estT(l)==0) || (rec.Alg.AlignedRec==2 && leff>=0) 
            EH.We=ones([size(yX,1) 1],'like',Residuals);
            for s=1:NSweeps;EH.We(mStSweeps==s)=We(s);end
        end
        if leff>0 && estT(l)==0
            if rec.Alg.WriteSnapshots && max(voxSiz)/min(voxSiz)<1.5;visSegment(flip(permute(dynInd(xRes,E.nF,3),[2 1 3]),2),[],2,1,[],[],sprintf('Corrected ($l=$%d)',l),strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.SuffOu,sprintf('_Di-l=%d',l)));end
            if parXT.computeCSRecon
                %save('fileNAMEHERE','yX','E','EH','P','CX','R','xRes','nX','tolType','estT','l','tol','-v7.3');   
                if nX~=0;xCS=CSsolver(yX,E,EH,P,CX,R,xRes,nX,tolType{estT(l)+1},tol(estT(l)+1),rec.Dyn.Debug);else xCS=xRes;end
                if rec.Alg.WriteSnapshots && max(voxSiz)/min(voxSiz)<1.5;visSegment(flip(permute(dynInd(xCS,E.nF,3),[2 1 3]),2),[],2,1,[],[],sprintf('Corrected regularized ($l=$%d)',l),strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.SuffOu,sprintf('_Re-l=%d',l)));end
            end
            EH.We=[];
        end

        if estT(l)               
            %LM SOLVER        
            %TRANSFORM PARAMETERS PER SEGMENT/HARMONICS/BLOCK SIZES            
            E.Tr=rec.T;E.Fs=FSt;
            E.NMs=size(E.Tr,5);E.NSe=length(E.Fs{1});
            if parXT.exploreMemory;E.cT(:)=1;end%To test flow without running the main methods       

            %EXTRACT RELEVANT ROI AND COILS AND CALL THE MOTION SOLVER                       
            nct=1:eivaS(2);
            E.Sf=dynInd(SRes,nct,4);
            ySf=dynInd(ySt,nct,4);
            E.dS=[E.NMs eivaS(2)];
            if gpuIn;E.Sf=gpuArray(E.Sf);end
            if parXT.fractionOrder~=0;E.Fd=GRes;end
            E=LMsolver(ySf,E,xRes,CT);
            rec.T=E.Tr;
            THist=cat(1,THist,gather(rec.T));
        else
            if l==L && parXT.saveFinal
                recW=rec;recW.Alg.SaveRaw=3;recW.Dyn.Batch=1;
                recW.E=E;recW.EH=EH;recW.P=P;recW.CX=CX;recW.R=R;recW.xRes=xRes;recW.yX=yX;recW.fullResid=fullResid;
                recW.timeState=timeState;recW.We=We;recW.Residuals=Residuals;recW.Sweeps=Sweeps;recW.iSt=iSt{l};
                recW.S=[];recW.y=[];
                writeRaw(recW);recW=[];
            end
            break;
        end
        E.Sf=[];

        %CHECK CONVERGENCE
        if all(E.cT(:))%Weak convergence
            %if l==L;estT(l)=0;else break;end%If last level we perform another estimation for x           
            estT(l)=0;            
            if rec.Alg.WriteSnapshots;visMotion(rec,voxSiz,timeState,2,strcat(folderSnapshots,filesep,'Motion'),strcat(fileName,rec.Plan.SuffOu,sprintf('_l=%d',l)));end                
        end
    end

    %RECONSTRUCTION TO BE USED TO INITIALIZE THE NEXT RESOLUTION LEVEL
    x=resampling(xRes,NX);
    if gpuIn;rec.M=gpuArray(rec.M);end
    if discMa;x=rec.M.*x;end
    
    if leff>0
        rec.d=dynInd(rec.d,leff,4,gather(x));
        if parXT.computeCSRecon
            if discMa;xCS=rec.M.*xCS;end
            rec.r=dynInd(rec.r,leff,4,gather(xCS));
        end
    end
    if ~isempty(En{l})
        EnSort=En{l};EnSort(iSt{l},:)=EnSort;
        ReSort=Re{l};ReSort(iSt{l},:)=ReSort;  
        mStSort=E.mSt;mStSort(iSt{l})=mStSort;
        rec.Par.Mine.DISORDER.Residuals{l}=gather(reshape(ReSort,[NYres([1:2 5]) size(ReSort,2)]));%%%PLAY WITH IST SO THAT THIS CAN BE MADE COMPARABLE AND ALSO LINKED TO MOTIONS
        rec.Par.Mine.DISORDER.Energy{l}=gather(reshape(EnSort,[NYres([1:2 5]) size(EnSort,2)]));%%%PLAY WITH IST SO THAT THIS CAN BE MADE COMPARABLE AND ALSO LINKED TO MOTIONS
        rec.Par.Mine.DISORDER.Motions{l}=THist;
        rec.Par.Mine.DISORDER.States{l}=gather(reshape(mStSort,NYres([1:2 5])));
    end

    if parXT.writeInter && leff>0 || l==L          
        E=[];EH=[];

        drec=rec.d;
        if parXT.computeCSRecon;rrec=rec.r;end
        ROIrec=rec.Enc.ROI;          
        for n=typ2RecI';datTyp=rec.Plan.Types{n};
            if gpu;rec.(datTyp)=gpuArray(rec.(datTyp));end
            rec.(datTyp)=dynInd(rec.(datTyp),1:leff,4);
            rec.(datTyp)=flip(rec.(datTyp),4);
        end
        if l==L && leff~=1;typ2RecL=setdiff(typ2Rec,[5;12]);elseif l==L && leff==1;typ2RecL=setdiff(typ2Rec,5);else typ2RecL=typ2RecI;end
        %FINAL ADJUSTMENTS
        for n=typ2RecL';datTyp=rec.Plan.Types{n};
            rec.(datTyp)=permute(rec.(datTyp),perm);
            rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,0,1);
        end
        for n=typ2RecI';datTyp=rec.Plan.Types{n};
            if rec.Alg.MargosianFilter;rec.(datTyp)=margosianFilter(rec.(datTyp),rec.Enc);end
            rec.(datTyp)=removeOverencoding(rec.(datTyp),rec.Alg.OverDec);%REMOVE OVERDECODING  
        end
        rec.Enc=rmfield(rec.Enc,'ROI');
        rec.Dyn.Typ2Wri(17)=1;%To write the motion estimates to file
        typ2RecI=setdiff(typ2RecI,12);
        writeData(rec);
        if l~=L
            rec.Enc.ROI=ROIrec;            
            rec.d=drec;
            if parXT.computeCSRecon;rec.r=rrec;end
        end
    end
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time till level %d: %.3f s\n',l,tend);end
    if l==L;break;end%Early termination if no outliered states detected
end
