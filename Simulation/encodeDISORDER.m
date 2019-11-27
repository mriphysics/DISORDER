function x=encodeDISORDER(x,S,etDir,etRefDir,A,B,y)

%ENCODEDISORDER   Forward model of the reconstruction in the presence of
%motion
%   X=ENCODEDISORDER(X,{S},{ETDIR},{ETREFDIR},{A},{B},{Y})  
%   * X is the image
%   * {S} is the coil array sensitivity map
%   * {ETDIR} is a cell array of motion factors (computed using
%   precomputeFactorsSincRigidTransform)
%   * {ETREFDIR} is a cell array of reference geometry factors (computed 
%   using precomputeFactorsSincRigidTransform)
%   * {A} contains the sampling scheme or the size of the Fourier space
%   * {B} contains the multiband sampling scheme
%   * {Y} contains the samples, for computing the residuals
%   ** X is the encoded data
%

if nargin<2;S=[];end
if nargin<3;etDir=[];end
if nargin<4;etRefDir=[];end
if nargin<5;A=[];end
if nargin<6 || isempty(B);B=cell(1,2);end
if nargin<7;y=[];end
gpu=isa(x,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end

NX=size(x);
if numel(A)==2;NY=A;
elseif ~isempty(A);NY=size(A);
else NY=NX;
end
NY(end+1:8)=1;
if ~isempty(y);NY(8)=size(y,8);end

if ~isempty(etDir);x=sincRigidTransform(x,etDir,1);elseif size(x,8)==1;x=repmat(x,[1 1 1 1 1 1 1 NY(8)]);end%T     
x=bsxfun(@times,x,S);%S
if ~isempty(etRefDir);x=sincRigidTransform(x,etRefDir,1);elseif size(x,8)==1;x=repmat(x,[1 1 1 1 1 1 1 NY(8)]);end%T

for m=1:2
    if ~isempty(B{m})
        if size(x,5+m)==1
            NR=ones(1,5+m);NR(5+m)=size(B{m},5+m);
            x=repmat(x,NR);
        end
        x=indDim(x,B{m},m);
    end  
    x=fold(x,m,size(x,m),NY(m));%U    
end
for m=1:2;x=fftGPU(x,m);end%F
if ~isempty(y);x=bsxfun(@minus,x,y);end
if ~isempty(A) && numel(A)~=2;x=bsxfun(@times,x,A);end%A