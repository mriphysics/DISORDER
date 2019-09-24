function x=decodeDISORDER(x,SH,etInv,etFovInv,A,B,M,P)

%DECODEDISORDER   Forward model of the reconstruction in the presence of
%motion
%   X=DECODEDISORDER(X,{SH},{ETINV},{ETFOVINV},{A},{B},{M},{P})  
%   * X is the image
%   * {SH} are the conjugated sensitivities
%   * {ETINV} is a cell array of motion factors (computed using
%   precomputeFactorsSincRigidTransform)
%   * {ETFOVINV} is a cell array of reference geometry factors (computed 
%   using precomputeFactorsSincRigidTransform)
%   * {A} contains the sampling scheme or the size of the Fourier space
%   * {B} contains the MB sampling scheme
%   * {M} applies a mask
%   * {P} applies a preconditioner
%   ** X is the encoded data
%

if nargin<2;SH=[];end
if nargin<3;etInv=[];end
if nargin<4;etFovInv=[];end
if nargin<5;A=[];end
if nargin<6 || isempty(B);B=cell(1,2);end
if nargin<7;M=[];end
if nargin<8;P=[];end

gpu=isa(x,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

NY=size(x);
if ~isempty(SH);NX=size(SH);else NX=NY;end

if ~isempty(A) && numel(A)~=2;x=bsxfun(@times,x,A);end%A^H
for m=1:2;x=ifftGPU(x,m,gpuF);end%F^H
for m=1:2    
    NR=ones(1,2);
    if ~isempty(B{m});NR(m)=size(B{m},5+m);end
    x=ifold(x,m,NX(m)/NR(m),NY(m));%U^H   
    if ~isempty(B{m})
        y=repmat(x,NR);y(:)=0;        
        y=indDim(y,B{m},m,x);
        x=y;
    end
end
if ~isempty(etFovInv);x=sincRigidTransform(x,etFovInv,0,[],[],0);end%T^H
x=sum(bsxfun(@times,x,SH),4);%S^H
if ~isempty(etInv);x=sincRigidTransform(x,etInv,0);else x=multDimSum(x,5:8);end%T^H
if ~isempty(M);x=bsxfun(@times,x,M);end%M
if ~isempty(P);x=bsxfun(@times,x,P);end%P