function x=ifctGPU(x,m,gpu,FH)

%IFCTGPU   Configurable GPU-based IDCT computation
%   X=IFCTGPU(X,M,GPU,{F})
%   * X is the array on which to apply the IDCT (Inverse Discrete Cosine
%   Transform)
%   * M is the direction along which to apply the IDCT
%   * GPU is a flag that determines whether to use cpu/gpu with built-in
%   matlab functions (0), gpu (1,3) or matrix-based gpu computation (2,4). 
%   These alternatives are introduced due to bugs in gpu computations in 
%   matlab 2015a for our architectures (roughly speaking these affected the 
%   gpu computation of the ifft for certain ---usually prime numbers--- 
%   sizes of the array when the computation was not performed along the 
%   first dimension), so users should better test the matlab gpu ifft 
%   behaviour in their systems with different array sizes and along 
%   different dimensions
%   * {F} is a IDCT matrix (or any other square matrix) provided by the 
%   user. If not provided and required, it is obtained by the function
%   * X is the IDCT-transformed array
%

if nargin<4;FH=[];end
if strcmp(version('-release'),'2015a') && gpu==0;gpu=1;end

ism=1;
if gpu==0
    x=ifct(x,[],m);
elseif gpu==3 || gpu==1
    if m==1
        if size(x,1)~=1;x=ifct(x);end
    else
        perm=1:ndims(x);perm(1)=m;perm(m)=1;
        x=permute(x,perm);
        if size(x,1)~=1;x=ifct(x);end
        x=permute(x,perm);
    end
elseif gpu==4 || gpu==2
    N=size(x,m);
    if N~=1
        if isempty(FH)
            if ~isaUnderlying(x,'double');FH=gpuArray(single(dctmtx(N)'));else FH=gpuArray(dctmtx(N)');end
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
end