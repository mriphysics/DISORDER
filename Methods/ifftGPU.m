function x=ifftGPU(x,m,FH,NX)

%IFFTGPU   Configurable GPU-based IDFT computation
%   X=IFFTGPU(X,M,{F},{NY})
%   * X is the array on which to apply the IDFT (Inverse Discrete Fourier
%   Transform)
%   * M is the direction along which to apply the IDFT
%   * {F} is a IDFT matrix (or any other square matrix) provided by the 
%   user. If not provided and required, it is obtained by the function
%   * {NX} are the desired spatial dimensions ---spatial truncation (NY>NX)
%   or unfolding (NY<NX)---
%   * X is the IDFT-transformed array
%

if nargin<3;FH=[];end
NY=size(x,m);
if nargin<4 || isempty(NX);NX=NY;end

gpu=isa(x,'gpuArray');
BlSz=1e3;
ND=ndims(x);

if NY~=1
   if isempty(FH);[~,FH]=build1DFTM(NY,0,gpu);end
   if NX<NY;FH=resampling(FH,[NX NY],2);end

    if (gpu && isaUnderlying(x,'double')) || isa(x,'double');FH=double(FH);end
    S=size(x);S(end+1:max(ND+1,m+1))=1;
    if m~=1;x=reshape(x,[prod(S(1:m-1)) S(m) prod(S(m+1:ND))]);else x=x(:,:);end
    if m==1
        x=FH*x;
    elseif m~=ND
        if numel(x)>1e8
            N3=size(x,3);N1=size(x,1);
            xo=resampling(x,[N1 size(FH,1)],2);
            for s=1:BlSz:N3;vS=s:min(s+BlSz-1,N3);xo(:,:,vS)=matfun(@mtimes,x(:,:,vS),FH.');end
            x=xo;xo=[];
        else
            x=matfun(@mtimes,x,FH.');
        end
    else
        x=x*FH.';
    end
    if m==1;S(m)=size(x,1);else S(m)=size(x,2);end
    x=reshape(x,S);
end

if NX>NY;x=ifold(x,m,NX,NY);end
