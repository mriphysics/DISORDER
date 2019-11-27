function x=fftGPU(x,m,F,NY)

%FFTGPU   Configurable GPU-based DFT computation
%   X=FFTGPU(X,M,{F},{NY})
%   * X is the array on which to apply the DFT (Discrete Fourier Transform)
%   * M is the direction along which to apply the DFT
%   * {F} is a DFT matrix (or any other square matrix) provided by the 
%   user. If not provided and required, it is obtained by the function
%   * {NY} are the desired spectral dimensions ---spatial zero pad (NY>NX) 
%   or folding (NY<NX)---
%   ** X is the DFT-transformed array
%

if nargin<3;F=[];end
NX=size(x,m);
if nargin<4 || isempty(NY);NY=NX;end

gpu=isa(x,'gpuArray');
BlSz=1e6;
ND=ndims(x);

if NY<NX
    x=fold(x,m,NX,NY);
    NX=NY;
end

if NY~=1
    if isempty(F);F=build1DFTM(NY,0,gpu);end
    if NX<NY;F=resampling(F,[NY NX],2);end
 
    if (gpu && isaUnderlying(x,'double')) || isa(x,'double');F=double(F);end
    S=size(x);S(end+1:max(ND+1,m+1))=1;
    if m~=1;x=reshape(x,[prod(S(1:m-1)) S(m) prod(S(m+1:ND))]);else x=x(:,:);end
    if m==1
        N2=size(x,2);
        xo=resampling(x,NY,2);
        for s=1:BlSz:N2;vS=s:min(s+BlSz-1,N2);xo(:,vS)=F*x(:,vS);end
        x=xo;xo=[];
    elseif m~=ND
        x=matfun(@mtimes,x,F.');
    else
        x=x*F.';
    end
    if m==1;S(m)=size(x,1);else S(m)=size(x,2);end            
    x=reshape(x,S);
end
