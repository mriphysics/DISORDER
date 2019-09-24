function x=ifftGPU(x,m,gpu,FH)

%IFFTGPU   Configurable GPU-based IDFT computation
%   X=IFFTGPU(X,M,GPU,{F})
%   * X is the array on which to apply the IDFT (Inverse Discrete Fourier
%   Transform)
%   * M is the direction along which to apply the IDFT
%   * GPU is a flag that determines whether to use cpu/gpu with built-in
%   matlab functions (0), gpu (1) or matrix-based gpu computation (2). 
%   These alternatives are introduced due to bugs in gpu computations in 
%   matlab 2015a for our architectures (roughly speaking these affected the 
%   gpu computation of the ifft for certain ---usually prime numbers--- 
%   sizes of the array when the computation was not performed along the 
%   first dimension), so users should better test the matlab gpu ifft 
%   behaviour in their systems with different array sizes and along 
%   different dimensions
%   * {F} is a IDFT matrix (or any other square matrix) provided by the 
%   user. If not provided and required, it is obtained by the function
%   * X is the IDFT-transformed array
%

if nargin<4;FH=[];end
if strcmp(version('-release'),'2015a') && gpu==0;gpu=1;end

BlSz=1e3;
if gpu==0 && isempty(FH)
    x=ifft(x,[],m);
elseif gpu==1 && isempty(FH)
     if m~=1
         perm=1:ndims(x);perm([1 m])=[m 1];
         x=permute(x,perm);
     end
     if size(x,1)~=1;x=ifft(x);end
     if m~=1;x=permute(x,perm);end
elseif gpu==2 || ~isempty(FH)
    N=size(x,m);   
    ND=ndims(x);
    if N~=1
        if isempty(FH);[~,FH]=build1DFTM(N,0,gpu);end
        if (gpu && isaUnderlying(x,'double')) || isa(x,'double');FH=double(FH);end
        S=size(x);S(end+1:max(ND+1,m+1))=1;
        if m~=1;x=reshape(x,[prod(S(1:m-1)) S(m) prod(S(m+1:ND))]);else x=x(:,:);end
        if m==1
            x=FH*x;
        elseif m~=ND
            if numel(x)>1e8
                N3=size(x,3);
                for s=1:BlSz:N3;vS=s:min(s+BlSz-1,N3);x(:,:,vS)=matfun(@mtimes,x(:,:,vS),FH.');end
            else
                x=matfun(@mtimes,x,FH.');
            end
        else
            x=x*FH.';
        end
        if m==1;S(m)=size(x,1);else S(m)=size(x,2);end
        x=reshape(x,S);
    end
end