function x=fctGPU(x,m,F)

%FCTGPU   Configurable GPU-based DCT computation
%   X=FCTGPU(X,M,{F})
%   * X is the array on which to apply the DCT (Discrete Cosine Transform)
%   * M is the direction along which to apply the DCT
%   * {F} is a FCT matrix (or any other square matrix) provided by the 
%   user. If not provided and required, it is obtained by the function
%   ** X is the FCT-transformed array
%

if nargin<3;F=[];end
gpu=isa(x,'gpuArray');

ism=1;
N=size(x,m);
if N~=1
    if isempty(F)
        F=dctmtx(N);
        if (gpu && ~isaUnderlying(x,'double')) || (~gpu && ~isa(x,'double'));F=single(F);end
        if gpu;F=gpuArray(F);end
    end
    if m~=1
       perm=1:ndims(x);perm([1 m])=[m 1];
       x=permute(x,perm);
    end
    if ~ismatrix(x)
        S=size(x);S(end+1:2)=1;
        x=reshape(x,[S(1) prod(S(2:end))]);
        ism=0;
    end
    x=F*x;
    if ~ism;x=reshape(x,S);end
    if m~=1;x=permute(x,perm);end
end
