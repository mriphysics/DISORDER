function x=fftGPU(x,m,gpu,F)

%FFTGPU   Configurable GPU-based DFT computation
%   X=FFTGPU(X,M,GPU,{F})
%   * X is the array on which to apply the DFT (Discrete Fourier Transform)
%   * M is the direction along which to apply the DFT
%   * GPU is a flag that determines whether to use cpu/gpu with built-in
%   matlab functions (0), gpu (1) or matrix-based gpu computation (2). 
%   These alternatives are introduced due to bugs in gpu computations in 
%   matlab 2015a for our architectures (roughly speaking these affected the 
%   gpu computation of the fft for certain ---usually prime numbers--- 
%   sizes of the array when the computation was not performed along the 
%   first dimension), so users should better test the matlab gpu fft 
%   behaviour in their systems with different array sizes and along 
%   different dimensions
%   * {F} is a DFT matrix (or any other square matrix) provided by the 
%   user. If not provided and required, it is obtained by the function
%   ** X is the DFT-transformed array
%

if nargin<4;F=[];end
if strcmp(version('-release'),'2015a') && gpu==0;gpu=1;end

BlSz=1e6;
if gpu==0 && isempty(F)
    x=fft(x,[],m);
elseif gpu==1 && isempty(F)
     if m~=1
         perm=1:ndims(x);perm([1 m])=[m 1];
         x=permute(x,perm);
     end
     if size(x,1)~=1;x=fft(x);end
     if m~=1;x=permute(x,perm);end
elseif gpu==2 || ~isempty(F)
    N=size(x,m);
    ND=ndims(x);        
    if N~=1
        if isempty(F);F=build1DFTM(N,0,gpu);end
        if (gpu && isaUnderlying(x,'double')) || isa(x,'double');F=double(F);end
        S=size(x);S(end+1:max(ND+1,m+1))=1;
        if m~=1;x=reshape(x,[prod(S(1:m-1)) S(m) prod(S(m+1:ND))]);else x=x(:,:);end
        if m==1
            N2=size(x,2);            
            for s=1:BlSz:N2;vS=s:min(s+BlSz-1,N2);x(:,vS)=F*x(:,vS);end
        elseif m~=ND
            x=matfun(@mtimes,x,F.');
        else
            x=x*F.';
        end
        if m==1;S(m)=size(x,1);else S(m)=size(x,2);end               
        x=reshape(x,S);
    end
end