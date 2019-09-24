function [et,etg,eth]=precomputeFactorsSincRigidTransform(kGrid,rkGrid,T,di,cg,gr,sh,cGrid,kkGrid)

%PRECOMPUTEFACTORSSINCRIGIDTRANSFORM precomputes the k-space phase multiplicative factors 
%required to apply a rigid transform based on sinc interpolation
%   [ET,ETG]=PRECOMPUTEFACTORSSINCRIGIDTRANSFORM(KGRID,RKGRID,T,{DI},{CG},{GR}) 
%   * KGRID is a grid of points in the spectral domain
%   * RKGRID is a grid of points in the spatial-spectral domain
%   * T are the parameters of the transform
%   * {DI} is a flag to indicate whether to perform direct or inverse 
%   transform (defaults to 1)
%   * {CG} is a flag to indicate whether to calculate the derivative terms 
%   (defaults to 0)
%   * {GR} is a flag to indicate whether this factors are for the grouped 
%   transform (defaults to 0)
%   * {SH} is a flag to indicate whether the grids have already been 
%   ifftshifted (defaults to 0)
%   * {CGRID} is the center of the grid
%   * {KKGRID} is a grid of points in the spectral-spectral domain
%   ** ET are the parameters to apply the transform
%   ** ETG are the factors to apply the derivative of the transform 
%   ** ETH are the factors to apply the Hessian of the transform 
%

if nargin<4 || isempty(di);di=1;end
if nargin<5 || isempty(cg);cg=0;end
if nargin<6 || isempty(gr);gr=0;end
if nargin<7 || isempty(sh);sh=0;end

N=zeros(1,length(kGrid));
for n=1:length(kGrid);N(n)=numel(kGrid{n});end
if nargin<8 || isempty(cGrid);cGrid=(N/2)+1;end
if (nargin<9 || isempty(kkGrid)) && cg==2
    fact=[1 2 3 1 1 2;
          1 2 3 2 3 3];
    kkGrid=cell(1,6);
    for m=1:6;kkGrid{m}=bsxfun(@times,kGrid{fact(1,m)},kGrid{fact(2,m)});end
end

per=[1 3 2;
     2 1 3];

%COMPUTATION OF LINEAR PHASE TERMS
et=cell(1,5);
ndT=ndims(T);
t=((-1)^di)*1i*dynInd(T,1:3,ndT);
theta=wrapToPi(dynInd(T,4:6,ndT));
gpu=isa(kGrid{1},'gpuArray');
for n=1:3
    et{5}{n}=dynInd(theta,n,ndT);et{5}{n}(:)=0;
    et{5}{n}(abs(dynInd(theta,n,ndT))>pi/2)=1;%These need to be flipped
    if gpu;et{5}{n}=gpuArray(et{5}{n});end
    for s=1:2;et{4}{n}{per(s,n)}=cGrid(per(s,n));end%Centers of the grid for flipping
end
theta=wrapToPiHalf(theta);

tantheta2=theta/2;
tantheta2=tan(tantheta2);
tantheta2j=((-1)^(di-1))*1i*tantheta2;
sintheta=sin(theta);
sintheta=((-1)^di)*1i*sintheta;
if cg>0
    if cg==2;tantheta=tan(theta);end
    tanthetacuad=tantheta2.*tantheta2;
    tanthetacuad=(1+tanthetacuad)/2;
    costheta=cos(theta);    
end       
theta=[];

%ASSIGNMENT 
if ~sh
    for m=1:3
        for n=1:2;rkGrid{n}{m}=ifftshift(rkGrid{n}{m},per(n,m));end
    end
end
 
et{2}=cell(1,3);
et{3}=cell(1,3);
for m=1:3    
    et{2}{m}=exp(bsxfun(@times,dynInd(tantheta2j,m,ndT),rkGrid{1}{m}));%Tan exponential
    et{3}{m}=exp(bsxfun(@times,dynInd(sintheta,m,ndT),rkGrid{2}{m}));%Sin exponential
end
tantheta2j=[];sintheta=[];

if cg>0
    etg=cell(1,3);
    etg{2}=cell(1,3);
    etg{3}=cell(1,3);
    for m=1:3     
        etg{2}{m}=1i*bsxfun(@times,dynInd(tanthetacuad,m,ndT),rkGrid{1}{m});%Tan derivative
        etg{3}{m}=-1i*bsxfun(@times,dynInd(costheta,m,ndT),rkGrid{2}{m});%Sin derivative
        if cg==2%Hessian terms
            eth{2}{m}=bsxfun(@plus,dynInd(tantheta2,m,ndT),etg{2}{m});
            eth{3}{m}=bsxfun(@plus,-dynInd(tantheta,m,ndT),etg{3}{m});
        end
    end
    tanthetacuad=[];costheta=[];tantheta2=[];tantheta=[];
    for m=2:3        
        for n=1:3
            etg{m}{n}=etg{m}{n}.*et{m}{n};
            if cg==2;eth{m}{n}=eth{m}{n}.*etg{m}{n};end
        end
    end
end

if ~sh
    for m=1:3
        kGrid{m}=ifftshift(kGrid{m},m);        
    end
end

if ~gr
    et{1}=exp(bsxfun(@plus,bsxfun(@plus,bsxfun(@times,dynInd(t,1,ndT),kGrid{1}),bsxfun(@times,dynInd(t,2,ndT),kGrid{2})),bsxfun(@times,dynInd(t,3,ndT),kGrid{3})));
else
    et{1}=cell(1,3);
    for m=1:3;et{1}{m}=exp(bsxfun(@times,dynInd(t,m,ndT),kGrid{m}));end
end
t=[];

if cg>0 
    etg{1}=cell(1,3);
    for m=1:3       
        if ~gr;etg{1}{m}=bsxfun(@times,et{1},-1i*kGrid{m});
        else etg{1}{m}=-1i*kGrid{m};
        end
    end
    if cg==2
        eth{1}=cell(1,6);
        for m=1:6
            if ~gr;eth{1}{m}=bsxfun(@times,-kkGrid{m},et{1});
            else eth{1}{m}=-kkGrid{m};end
        end
    end
end
