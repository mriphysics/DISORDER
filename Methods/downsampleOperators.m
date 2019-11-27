function [xRes,SRes,yRes,ARes,MRes,rGridRes,kGridRes,rkGridRes,F,FH,fail]=downsampleOperators(l,x,S,y,A,M,gibbsRing,fo,gpuIn)

%DOWNSAMPLEOPERATORS   Spatially downsamples the aligned reconstruction operators
%   [XRES,SRES,yRES,ARES,MRES,RGRIDRES,KGRIDRES,RKGRIDRES]=DOWNSAMPLEOPERATORS(L,X,S,Y,A,M,RGRID,KGRID,RKGRID,{GIBBSRING},{FO},{GPUIN})
%   * L is the number of 2x subdivisions
%   * X is the image
%   * S are the sensitivity maps
%   * Y are the measurements
%   * A is the sampling operator
%   * M is the mask
%   * {GIBBSRING} indicates whether to use Gibbs ringing filtering
%   * {FO} indicates if y is in Fourier domain
%   * {GPUIN} serves to indicate the gpu nature of the transform and Fourier operators
%   ** X is the resampled image
%   ** S are the resampled sensitivity maps
%   ** Y are the resampled measurements
%   ** A is the resampled sampling operator
%   ** M is the resampled mask
%   ** RGRID is the resampled spatial grid
%   ** KGRID is the resampled spectral grid
%   ** RKGRID is the resampled spatial-spectral grid
%   ** F are the resampled Fourier operators
%   ** FH are the resampled back-Fourier operators
%   ** FAIL indicates that the reconstruction has to be aborted due to insuficient sampling
%

if nargin<7 || isempty(gibbsRing);gibbsRing=0;end
if nargin<8 || isempty(fo);fo=1;end

gpu=isa(x,'gpuArray');
if nargin<9 || isempty(gpuIn);gpuIn=gpu;end

NX=size(x);NY=size(y);NX(end+1:3)=1;NY(end+1:12)=1;
NXRes=ceil(NX(1:3)./(2.^l));NYRes=ceil(NY(1:3)./(2.^l));

fail=0;
ARes=resampling(A,NYRes(1:2),1);
if any(multDimSum(ARes,1:2)==0)%For sequential we can't accelerate along the first phase encode
    fail=1;
    fprintf('Sampling does not seem to be distributed, suboptimal performance is to be expected\n');
    NXRes(1)=NX(1);NYRes(1)=NY(1);
    ARes=resampling(A,NYRes(1:2),1);        
end

xRes=resampling(x,NXRes);
NS=size(S);NS(end+1:4)=1;
SRes=zeros([NXRes(1:3) NS(4)],'like',S);
yRes=zeros([NYRes(1:3) NS(4) NY(5:12)],'like',y);
for s=1:NS(4)
    SRes=dynInd(SRes,s,4,resampling(dynInd(S,s,4),NXRes));
    yRes=dynInd(yRes,s,4,resampling(dynInd(y,s,4),NYRes,fo));
end

if gibbsRing~=0 && gibbsRing<=1;HY=buildFilter(NYRes(1:3),'tukeyIso',[],gpuIn,gibbsRing);elseif l~=0;HY=buildFilter(NYRes,'CubicBSpline',[],gpuIn);else HY=[];end
if ~isempty(HY)
    if ~fo;yRes=filtering(yRes,HY);else yRes=bsxfun(@times,yRes,HY);end
end

MRes=resampling(M,NXRes);
if all(ismember(M(:),[0 1]));MRes=single(MRes>0.5);else MRes=abs(MRes);end%Discrete mask    

centRes=ceil((NXRes+1)/2);
[rGridRes,kGridRes,rkGridRes]=generateTransformGrids(NX,gpuIn,NXRes,centRes,1);

[F,FH]=buildStandardDFTM(NXRes(1:3),0,gpuIn);
