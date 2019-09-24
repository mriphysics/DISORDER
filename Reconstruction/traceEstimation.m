function [tr,trF]=traceEstimation(x,E,EH,C,R,tol,BlSz)

%TRACEESTIMATION   Estimates the traces of operators involved in image
%reconstruction by the normalized Rayleigh-quotient trace estimator in [1] 
%H Avron, S Toledo, "Randomized algorithms for estimating the trace of an 
%implicit symmetric positive semi-definite matrix", J ACM, 58(2):8:1-13, 
%Apr 2011.
%   TR=TRACEESTIMATION(X,E,C,R,TOL,BLSZ)
%   * X is a variable from which to infer the matrix dimensions
%   * {E} is an encoding structure
%   * {EH} is a decoding structure
%   * {C} is a constrain structure
%   * {R} is a regularizer structure
%   * {TOL} is the ratio of the variance of the estimate with regard to the
%   estimated value (for convergence). It defaults to 1e-2
%   * {BLSZ} is the number of realizations of parallel computations. It 
%   defaults to 1
%   ** TR is the estimated trace
%   ** TRF are the set of estimates for the trace
%

if nargin<2;E=[];end
if nargin<3;EH=[];end
if nargin<4;C=[];end
if nargin<5;R=[];end
if nargin<6 || isempty(tol);tol=1e-2;end
if nargin<7 || isempty(BlSz);BlSz=1;end

gpu=isa(x,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end
x=repmat(real(x),[1 1 1 1 1 BlSz]);
NX=size(x);NX(end+1:6)=1;NX(6)=BlSz;
N=prod(NX)/BlSz;
conv=0;

if ~isempty(E)
    if ~isfield(E,'bS')%Block sizes by default
        E.bS=[size(x,5) size(x,4)];E.dS=[size(x,5) size(x,4)];E.oS=[1 1];
        if isfield(E,'Sf');E.bS(2)=size(E.Sf,4);E.dS(2)=size(E.Sf,4);end
        if isfield(E,'Tr');E.bS(1)=E.NSe;E.dS(1)=E.NSe;end
    end
end
if ~isempty(EH)
    if isfield(EH,'We')
        iSt=cat(1,E.iSt{:});
        EH.We(iSt)=EH.We;
        NWe=numel(EH.We);
        EH.We=reshape(EH.We,[E.Uf{1}.NY E.Uf{2}.NY 1 1 NWe/(E.Uf{1}.NY*E.Uf{2}.NY)]);
    end
end     

NitTot=0;
trF=[];
while ~conv
    x=plugNoise(x,1);
    x=bsxfun(@rdivide,x,sqrt(normm(x,[],1:5)/N));
    trC=zeros([1 1 1 1 1 BlSz],'like',x);
    %CONSTRAIN
    if ~isempty(C)
        if isfield(C,'Ma');x=bsxfun(@times,x,C.Ma);end
    end
    if ~isempty(E)
        %RESAMPLING
        if isfield(E,'NXAcq');x=resampling(x,E.NXAcq);end            
        %COIL PROFILES
        for b=E.oS(2):E.bS(2):E.dS(2);vB=b:min(b+E.bS(2)-1,E.dS(2));
            xS=x;
            if isfield(E,'Sf')
                Saux=dynInd(E.Sf,vB,4);
                if isa(xS,'gpuArray');Saux=gpuArray(Saux);end
                xS=bsxfun(@times,xS,Saux);
            end%Sensitivities
            %SENSE
            if isfield(E,'Uf')%Sense folding (first two dimensions)
                for n=1:2
                    if ~isempty(E.Uf{n});xS=fold(xS,n,E.Uf{n}.NX,E.Uf{n}.NY,E.UAf{n});end
                end
            end
            if isfield(EH,'We')
                for n=1:2;xS=fftGPU(xS,n,gpuF)/sqrt(size(xS,n));end
                xS=bsxfun(@times,xS,sqrt(EH.We));
            end            
            trC=trC+normm(xS,[],1:5);
        end
    end
    if ~isempty(R)
        %SHEARLET
        ND=3;%Dimensions
        if isfield(R,'Sh')            
            xF=x;
            for m=1:ND;xF=fftGPU(xF,m,gpuF);end
            NW=size(R.Sh.We,4);
            for w=1:NW
                sh=R.Sh.sH.S(:,:,:,w);
                we=R.Sh.We(:,:,:,w);
                if gpu;sh=gpuArray(sh);we=gpuArray(we);end
                xS=bsxfun(@times,xF,conj(sh));
                for m=1:ND;xS=ifftGPU(xS,m,gpuF);end
                xS=bsxfun(@times,xS,sqrt(we));
                trC=trC+normm(xS,[],1:4);
            end
        end        
    end    
    trC=gather(trC);
    trF=cat(6,trF,trC);
    NitTot=NitTot+BlSz;
    trMea=mean(trF,6);
    trStd=std(trF,[],6)/sqrt(NitTot);
    if NitTot>=2 && (trStd/(trMea+eps))<tol
        tr=trMea;
        conv=1;
        fprintf('Number of samples for trace estimation: %d\n',NitTot);
    end    
end
