function x=ifctGPU(x,m,FH)

%IFCTGPU   Configurable GPU-based IDCT computation
%   X=IFCTGPU(X,M,{F})
%   * X is the array on which to apply the IDCT (Inverse Discrete Cosine
%   Transform)
%   * M is the direction along which to apply the IDCT
%   * {F} is a IDCT matrix (or any other square matrix) provided by the 
%   user. If not provided and required, it is obtained by the function
%   * X is the IDCT-transformed array
%

if nargin<3;FH=[];end
gpu=isa(x,'gpuArray');

ism=1;
N=size(x,m);
if N~=1
    if isempty(FH)
        FH=dctmtx(N)';
        if (gpu && ~isaUnderlying(x,'double')) || (~gpu && ~isa(x,'double'));FH=single(FH);end
        if gpu;FH=gpuArray(FH);end
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
    x=FH*x;
    if ~ism;x=reshape(x,S);end
    if m~=1;x=permute(x,perm);end
end
