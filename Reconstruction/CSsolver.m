function [x,lambda]=CSsolver(y,E,EH,A,C,R,x,nIt,tolType,tol,deb,nOu)

%CSSOLVER   Performs a CG-based CS pseudoinverse reconstruction
%   [X,EN]=CSSOLVER(Y,E,EH,{A},{C},{R},{X},{NIT},{TOLTYPE},{TOL},{DEB},{NOU})
%   * Y is the array of measures
%   * E is an encoding structure
%   * EH is a decoding structure
%   * {A} is a preconditioner structure
%   * {C} is a constrain structure
%   * {R} is a regularizer structure
%   * {X} is an initial condition for the array to be reconstructed
%   * {NIT} is the maximum number of iterations
%   * {TOLTYPE} is the type of tolerance used for convergence
%   * {TOL} is the tolerance set for convergence
%   * {DEB} indicates whether to print information about convergence
%   * {NOU} indicates the number of outer iterations for CS reconstruction
%   ** X is the reconstruction
%   ** LAMBDA is the estimated lambda
%

if nargin<4;A=[];end
if nargin<5;C=[];end
if nargin<6;R=[];end
if nargin<8 || isempty(nIt);nIt=300;end%300;end
if nargin<9 || isempty(tolType);tolType='Energy';end%tolType='NormwiseBackward2Error'
if nargin<10 || isempty(tol);tol=1e-5;end%5e-3;end%1e-2;end
if nargin<11 || isempty(deb);deb=1;end
if nargin<12 || isempty(nOu);nOu=2;end

if deb>=1;tsta=tic;end

%INITIALIZATION
if isfield(EH,'Mbe') && ~isempty(EH.Mbe);EH.Mbe(:)=1;end
if nargin<7 || isempty(x)
    if isfield(EH,'We');We=EH.We;EH.We(:)=1;end
    x=CGsolver(y,E,EH,A,C,R,[],nIt,tolType,tol,[],[],deb);
    if isfield(EH,'We');EH.We=We;end
end
gpu=isa(x,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end

%RESAMPLING
NN=size(x);
E.NXAcq=NN;
EH.NXRec=NN;
x=resampling(x,EH.NXRec);
if ~isempty(C);C.Ma=single(resampling(C.Ma,EH.NXRec)>0.5);
else R.Ti.la=single(resampling(R.Ti.la,EH.NXRec));%THIS MAY BE PROBLEMATIC
end

%BUILD SHEARLET
epsil=1e-2;%To condition
R.Sh.p=1;%L1 norm
J=1;%1 level
R.Sh.sH=buildShearlet(NN,J,0,ones(1,J),'kos');
NW=size(R.Sh.sH.S,4);
if deb>=1;fprintf('Size shearlets:%s %d\n',sprintf(' %d',NN),NW);end
R.Sh.We=zeros([NN NW],'like',R.Sh.sH.S);
computeShearlet;
R.Sh.We=(R.Sh.We+epsil).^(R.Sh.p-2);

%ESTIMATE LAMBDA
lambda=zeros(1,nOu+1);
trFid=traceEstimation(x,E,EH,C,[],1e-2,4);
trReg=traceEstimation(x,[],[],C,R,1e-2,4);
R.Sh.la=(2/R.Sh.p)*trFid/(trReg+eps);
lambda(1)=R.Sh.la;
fprintf('It 0. Lambda: %.6g\n',R.Sh.la);

for n=1:nOu
    %PRECONDITIONER
    A.Se=normm(E.Sf,[],4);     
    if isempty(C);A.Se=A.Se+R.Ti.la;end
    A.Se=(A.Se+R.Sh.la*R.Sh.p/2).^(-1);
    A.Se=resampling(A.Se,EH.NXRec);
     
    x=CGsolver(y,E,EH,A,C,R,x,3,'None',tol);
    if n~=nOu
        computeShearlet;
        %COMPUTE SHEARLETS
        R.Sh.We=(R.Sh.We+epsil).^(R.Sh.p-2);
    
        %COMPUTE NEW LAMBDA
        trReg=traceEstimation(x,[],[],C,R,1e-2,4);
        R.Sh.la=(2/R.Sh.p)*trFid/(trReg+eps);
        lambda(n+1)=R.Sh.la;
        fprintf('It %d. Lambda: %.6g\n',n,R.Sh.la);
    end
end   
if deb>=1;tend=toc(tsta);fprintf('Time for final recon: %.3f s\n',tend);end

function computeShearlet
    xF=x;
    for m=1:3;xF=fftGPU(xF,m,gpuF);end
    for w=1:NW
        sh=R.Sh.sH.S(:,:,:,w);
        if gpu;sh=gpuArray(sh);end
        sh=bsxfun(@times,conj(sh),xF);
        for m=1:3;sh=ifftGPU(sh,m,gpuF);end
        R.Sh.We(:,:,:,w)=gather(abs(sh));
    end
end

end
