function rec=solveXT(rec)

%SOLVEXT   Performs an aligned reconstruction according to the manuscript
%"Motion corrected magnetic resonance imaging with DISORDER: Distributed 
%and Incoherent Sample Orders for Reconstruction Deblurring using Encoding 
%Redundancy", L Cordero-Grande, G Ferrazzi, RPAG Teixeira, AN Price, 
%J O’Muircheartaigh, and JV Hajnal, 2019.
%   REC=SOLVEXT(REC)
%   * REC is a reconstruction structure with input information
%   ** REC is a reconstruction structure with reconstructed data and
%   surrogate information
%

tsta=tic;

%SHORTCUTS FOR COMMONLY USED STRUCTURE FIELDS AND INITIALIZERS
voxSiz=rec.Enc.VoxelSize;%Acquired voxel size
parXT=rec.Alg;
gpu=single(gpuDeviceCount && ~blockGPU);if gpu;gpuF=2;else gpuF=0;end
on=cell(1,4);for n=1:4;on{n}=ones(1,n);end
typ2Rec=rec.Dyn.Typ2Rec;

%RECONSTRUCTION PLAN
[resPyr,L,estT,resIso]=pyramidPlan(voxSiz,parXT.resolMax,parXT.nwEnd,parXT.accel);

%PERMUTE PHASE ENCODES TO THE FIRST DIMENSIONS
perm=1:rec.Dyn.NDims;perm(1:3)=[3 2 1];
for n=typ2Rec'
    if ~ismember(n,[5 12]);datTyp=rec.Dyn.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
end
voxSiz=voxSiz(perm(1:3));parXT.apod=parXT.apod(perm(1:3));

%APODIZE + REARRANGE DATA + CORRECT FOR INTERLEAVED REPEATS
if gpu;rec.M=gpuArray(rec.M);end
if gpu;rec.y=gpuArray(rec.y);end
NY=size(rec.y);NY(end+1:rec.Dyn.NDims)=1;
rec.y=bsxfun(@times,rec.y,ifftshift(buildFilter(NY(1:3),'tukey',[10 10 1],gpu,parXT.apod)));%Apodize
[rec.y,NY]=resSub(rec.y,5:rec.Dyn.NDims);NY(end+1:5)=1;%Different averages along the 5th dimension
if gpu;rec.y=gather(rec.y);end

%ACQUISITION STRUCTURE
NEchos=rec.Enc.NEchos;
NProfs=numel(rec.Enc.K{2});
if NProfs<=1 || NProfs~=numel(rec.Enc.K{3})
    fprintf('SEPARABLE Y-Z TRAJECTORIES. ALIGNED RECONSTRUCTION IS NOT PERFORMED\n')
    if nargout<1;rec=[];else rec.Fail=1;end
    return;
end
if sum(cdfFilt(abs(diff(rec.Enc.K{2}(:))),'med'))<sum(cdfFilt(abs(diff(rec.Enc.K{3}(:))),'med'));fprintf('Slow direction is: %s\n',rec.Enc.MPS(4:5));else fprintf('Slow direction is: %s\n',rec.Enc.MPS(7:8));end
kTraj=zeros([NProfs 2],'single');
for n=1:2;kTraj(:,n)=rec.Enc.K{perm(n)}(:);end
if NEchos==1;NEchos=NProfs;end
NShots=NProfs/NEchos;
if mod(NShots,1)~=0
    fprintf('NUMBER OF SHOTS: %.2f, PROBABLY AN INTERRUPTED SCAN\n',NShots);
    if nargout<1;rec=[];else rec.Fail=1;end
    return;
end
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
%PLOT TRAJECTORIES
kRange=single(zeros(2,2));
for n=1:2;kRange(n,:)=rec.Enc.KRange{perm(n)};end
kShift=(floor((diff(kRange,1,2)+1)/2)+1)';
kSize=(diff(kRange,1,2)+1)';
kIndex=bsxfun(@plus,kTraj,kShift);
rec.DISORDER.kTraj=gather(kTraj);
rec.DISORDER.kRange=gather(kRange);
rec.DISORDER.kShift=gather(kShift);
rec.DISORDER.kSize=gather(kSize);
rec.DISORDER.kIndex=gather(kIndex);

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
    NSweeps=round(NProfs/NSamples);
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
for n=1:length(kIndex)
    for s=1:NY(5)
        if timeIndex(kIndex(n,1),kIndex(n,2),s)==0
            timeIndex(kIndex(n,1),kIndex(n,2),s)=n;
            break;
        end
    end
end
for m=1:2;timeIndex=ifftshift(timeIndex,m);end
if isSteadyState
    timeSample=(0:NProfs-1)*rec.Enc.TR(1)/1000;
else
    TPerShot=rec.Enc.T/NShots;    
    timeSample=(sweepSample-1)*TPerShot;
    if strcmp(rec.Enc.Technique,'TSE') || strcmp(rec.Enc.Technique,'TIR')
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*(rec.Enc.TE(1)/(1000*NEchos/2));
    else
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*rec.Enc.TR(1)/1000;
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
x=zeros(NX,'single');if gpu;x=gpuArray(x);end
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
iSt=cell(1,L);
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
[E.UAf,EH.UAb]=buildFoldM(NX(1:3),NY(1:3),gpu,1);%Folding/Unfolding matrix (empty is quicker!)
R=[];%Regularizer
ncx=1:parXT.eivaS(1);

%MASKING
CX.Ma=rec.M;
if parXT.softMasking==2;CX.Ma(:)=1;end
discMa=all(ismember(CX.Ma(:),[0 1]));
if ~discMa%Soft masking
    M=buildFilter(2*[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX],'tukeyIso',1,gpu,1,1);
    Ti.la=abs(filtering(CX.Ma,M,1));M=[];CX=[];
    R.Ti.la=1./(abs(Ti.la).^2+0.01);
end

%SENSITIVITIES
if gpu;rec.S=gpuArray(rec.S);end
E.Sf=dynInd(rec.S,ncx,4);E.dS=[1 ncx(end)];
if gpu;rec.S=gather(rec.S);end

%PRECONDITION
if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
else P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
end

if gpu;rec.y=gpuArray(rec.y);end
yX=mean(dynInd(rec.y,ncx,4),5);
if gpu;rec.y=gather(rec.y);end
if parXT.gibbsRinging~=0 && parXT.gibbsRinging<=1;yX=filtering(yX,buildFilter(NY(1:3),'tukeyIso',[],gpu,parXT.gibbsRinging));end
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
end
E.Sf=[];yX=[];

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
    if leff<=1
        if ~isSteadyState;timeState=cat(3,timeStateLimitInf,timeStateLimitSup);else timeState=cat(3,timeStateLimitInf,timeState,timeStateLimitSup);end%NOTE THIS ONLY WORKS IF THERE IS NO INTERSHOT MOTION
    end
   
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

    %RELATIVE SPATIAL RESOLUTION CALCULATIONS AND BLOCK SIZES FOR GPU COMPUTATION AT THIS LEVEL
    resAni=single(on{3});
    [~,indMinRes]=min(voxSiz);        
    resAni(indMinRes)=voxSiz(indMinRes)/resPyr(l);
    indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
    resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes));   
    resAni=voxSiz./resAni;
    
    %ENCODING STRUCTURES AND BLOCK SIZES FOR GPU COMPUTATION AT THIS LEVEL
    if gpu;rec.S=gpuArray(rec.S);rec.y=gpuArray(rec.y);end
    if ~discMa
        [xRes,SRes,yRes,timeIndexRes,R.Ti.la,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,rec.S,rec.y,timeIndex,Ti.la,parXT.gibbsRinging,0,gpu);
        R.Ti.la=1./(abs(R.Ti.la).^2+0.01);
    else
        [xRes,SRes,yRes,timeIndexRes,CX.Ma,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,rec.S,rec.y,timeIndex,rec.M,parXT.gibbsRinging,0,gpu);
    end
    if gpu;rec.S=gather(rec.S);rec.y=gather(rec.y);end
    
    NXres=size(xRes);NXres(end+1:3)=1;NYres=size(yRes);NYres(end+1:5)=1;
    for n=1:2;yRes=fftGPU(yRes,n,gpuF)/sqrt(NYres(n));end
    
    %WEIGHTING FUNCTION FOR MOTION ESTIMATION (ROI)
    if strcmp(rec.Enc.FatShift,'F');E.nF=1:floor((1-parXT.redFOV)*NXres(3));elseif strcmp(rec.Enc.FatShift,'H');E.nF=floor(1+parXT.redFOV*NXres(3)):NXres(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NXres(3);end
    EH.Mbe=zeros([1 1 NYres(3)],'like',real(xRes));EH.Mbe(E.nF)=1;%To use the extracted region for energy computation
    
    %HARMONIC SAMPLING
    E.mSt=timeIndexRes;  
    E.mSt(timeIndexRes~=0)=stateSample(timeIndexRes(timeIndexRes~=0));

    if ~isfield(EH,'We');EH.We=[];end    
    [ySt,FSt,EH.We,E.iSt,E.nSt,E.mSt,nSa]=buildHarmonicSampling(yRes,EH.We,E.mSt,NStates);yRes=[];iSt{l}=cat(1,E.iSt{:});%Samples as cell arrays
    if any(nSa(1:NStates)==0);nowithin=1;continue;end
        
    if leff<=1;mStSweeps=E.mSt;end
    BlSz=max(sqrt(prod(voxSiz(1:2)./resAni(1:2))),1);%Block sizes baseline
    E.bS=round([2 4].^log2(BlSz));E.oS=on{2};%Block sizes for our GPU
    if rec.Dyn.Debug>=2;fprintf('Block sizes:%s\n',sprintf(' %d',E.bS));end
    E.Je=1;%To use joint encoding/decoding
    E.Uf=cell(1,3);for n=1:3;E.Uf{n}.NY=NYres(n);E.Uf{n}.NX=NXres(n);end%Folding
    EH.Ub=cell(1,3);for n=1:3;EH.Ub{n}.NY=NYres(n);EH.Ub{n}.NX=NXres(n);end%Unfolding
    [E.UAf,EH.UAb]=buildFoldM(NXres(1:3),NYres(1:3),gpu,1);%Folding/Unfolding matrix (empty is quicker!)
    
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
        if (n==1 || l<L) && any(subTT(:)==2);E.cT(subTT==1)=1;end%For within-shot in the first iteration we only estimate those coming from outliers
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
        ncx=1:parXT.eivaS(1+2*estT(l));
        E.Sf=dynInd(SRes,ncx,4);E.dS(2)=ncx(end);
        yX=dynInd(ySt,ncx,4);

        %PRECONDITION
        if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
        else P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
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
        yX=dynInd(ySt,ncx,4);
        Re{l}(:,n)=computeEnergy(yX,xRes,E,[],EH,[],[],[],[],1)/numel(yX);
        En{l}(:,n)=normm(yX,[],setdiff(1:numDims(yX),1))/numel(yX);
        if rec.Dyn.Debug>=2;fprintf('Residuals: %0.6g\n',sum(Re{l}(:,n)));end        
        
        %ENERGY PER STATES
        Residuals=Re{l}(:,n);Residuals(iSt{l})=Residuals;
        Residuals=reshape(Residuals,NYres([1:2 5]));
        mStSort=mStSweeps;mStSort(iSt{l})=mStSort;
        Sweeps=reshape(mStSort,NYres([1:2 5]));         
        WeO=zeros([1 NSweeps],'like',Residuals);
        for s=1:NSweeps;WeO(s)=median(Residuals(Sweeps==s));end                      
        Westd=1.1926*multDimMed(abs(bsxfun(@minus,WeO',WeO)),1:2);%From  Rousseeuw, Peter J.; Croux, Christophe (December 1993), "Alternatives to the Median Absolute Deviation", Journal of the American Statistical Association, American Statistical Association, 88 (424): 1273–1283, doi:10.2307/2291267, JSTOR 2291267
        Wemea=sort(WeO);       
        Wemea=Wemea(min(max(round(parXT.percRobustShot*NSweeps+0.5),1),NSweeps));
        Wemea=Wemea+Westd*sqrt(2)*erfcinv(1-2*parXT.percRobustShot);        
        We=WeO/Wemea;
        We=1./We;  
        outlWe=We<parXT.enerRobustShot;
        
        %CS-LIKE/ROBUST RECONSTRUCTION
        if leff>0 && estT(l)==0
            subT(:)=1;
            if parXT.alignedRec~=4;outlWeV=find(outlWe);else outlWeV=1:length(outlWe);end
            for o=1:length(outlWeV);subT(sweepToSample{outlWeV(o)})=2;end
                
            if all(subT(:)==1) || nowithin;L=l;end
            EH.We=ones([size(yX,1) 1],'like',Residuals);%ALL THIS COULD BE PROBLEMATIC IN CASE OF SEVERAL REPEATS
            if ismember(parXT.alignedRec,2:3)
                for s=1:NSweeps;EH.We(mStSweeps==s)=1-single(outlWe(s));end
            else
                for s=1:NSweeps;EH.We(mStSweeps==s)=We(s);end
            end
            if parXT.computeCSRecon
                if nX~=0;xCS=CSsolver(yX,E,EH,P,CX,R,xRes,nX,tolType{estT(l)+1},tol(estT(l)+1),rec.Dyn.Debug);else xCS=xRes;end
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
            nct=1:parXT.eivaS(2);
            E.Sf=dynInd(SRes,nct,4);
            ySf=dynInd(ySt,nct,4);
            E.dS=[E.NMs parXT.eivaS(2)];
            E=LMsolver(ySf,E,xRes,CT);
            rec.T=E.Tr;
            THist=cat(1,THist,gather(rec.T));
        else
            break;
        end
        E.Sf=[];
        
        %CHECK CONVERGENCE
        if all(E.cT(:));estT(l)=0;end%Weak convergence                         
    end
    
    %RECONSTRUCTION TO BE USED TO INITIALIZE THE NEXT RESOLUTION LEVEL
    x=resampling(xRes,NX);
    if discMa;x=rec.M.*x;end
    if leff>0
        rec.d=dynInd(rec.d,leff,4,gather(x));
        if parXT.computeCSRecon
            if discMa;xCS=rec.M.*xCS;end
            rec.r=dynInd(rec.r,leff,4,gather(xCS));
        end
    end
    EnSort=En{l};EnSort(iSt{l},:)=EnSort;
    ReSort=Re{l};ReSort(iSt{l},:)=ReSort;  
    mStSort=E.mSt;mStSort(iSt{l})=mStSort;
    rec.DISORDER.Residuals{l}=gather(reshape(ReSort,[NYres([1:2 5]) size(ReSort,2)]));%%%PLAY WITH IST SO THAT THIS CAN BE MADE COMPARABLE AND ALSO LINKED TO MOTIONS
    rec.DISORDER.Energy{l}=gather(reshape(EnSort,[NYres([1:2 5]) size(EnSort,2)]));%%%PLAY WITH IST SO THAT THIS CAN BE MADE COMPARABLE AND ALSO LINKED TO MOTIONS
    rec.DISORDER.Motions{l}=THist;
    rec.DISORDER.States{l}=gather(reshape(mStSort,NYres([1:2 5])));

    if parXT.writeInter && leff>0 || l==L          
        E=[];EH=[];

        drec=rec.d;
        if parXT.computeCSRecon;rrec=rec.r;end
        for n=typ2RecI';datTyp=rec.Dyn.Types{n};
            if gpu;rec.(datTyp)=gpuArray(rec.(datTyp));end
            rec.(datTyp)=dynInd(rec.(datTyp),1:leff,4);
            rec.(datTyp)=flip(rec.(datTyp),4);
        end
        if l==L && leff~=1;typ2RecL=setdiff(typ2Rec,[5;12]);elseif l==L && leff==1;typ2RecL=setdiff(typ2Rec,5);else typ2RecL=typ2RecI;end
        %FINAL ADJUSTMENTS
        for n=typ2RecL';datTyp=rec.Dyn.Types{n};
            rec.(datTyp)=permute(rec.(datTyp),perm);
        end
        for n=typ2RecI';datTyp=rec.Dyn.Types{n};
            if parXT.margosianFilter;rec.(datTyp)=margosianFilter(rec.(datTyp),rec.Enc);end
            rec.(datTyp)=removeOverencoding(rec.(datTyp),rec.Enc.OverDec);%REMOVE OVERDECODING  
        end
        rec.Dyn.Typ2Wri(17)=1;%To write the motion estimates to file
        typ2RecI=setdiff(typ2RecI,12);
        writeData(rec);
        if l~=L
            rec.d=drec;
            if parXT.computeCSRecon;rec.r=rrec;end
        end
    end
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time till level %d: %.3f s\n',l,tend);end
    if l==L;break;end%Early termination if no outliered states detected
end
if gpu;rec.M=gather(rec.M);end
if nargout<1;rec=[];end