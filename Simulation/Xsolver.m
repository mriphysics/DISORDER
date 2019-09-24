function [x,n,E]=Xsolver(x,y,S,M,T,TFov,A,B,kGrid,rkGrid,nX,toler,lambda)

%XSOLVER   Reconstructs an image using CG SENSE with multishot alignment
%   X=XSOLVER(X,Y,S,{M},{T},{A},{KGRID},{RKGRID},{NX},{TOLER},{LAMBDA}) 
%   computes the best X for a given T
%   * X is the image to be reconstructed
%   * Y is the measured data
%   * S is the coil-array sensitivity map
%   * {M} is a spatial mask to constrain the solution
%   * {T} are the transform parameters
%   * {TFOV} are the transform parameters to the acquisition geometry
%   * {A} is a sampling mask
%   * {B} is the multiband sampling pattern
%   * {KGRID} is a grid of points in the spectral domain
%   * {RKGRID} is a grid of points in the spatial-spectral domain
%   * {NX} is the number of iterations of the CG algorithm
%   * {TOLER} is the maximum update for convergence. It defaults to 0
%   * {LAMBDA} is a regularization parameter. It defaults to 0
%   ** X the reconstructed image
%   ** N is the number of iterations till convergence
%   ** E is the energy of the residuals
%

NY=size(y);
if nargin<4;M=[];end
if nargin<5;T=[];end
if nargin<6;TFov=[];end
if nargin<7;A=NY(1:2);end
if nargin<8 || isempty(B);B=cell(1,2);end
if nargin<9;kGrid=[];end
if nargin<10;rkGrid=[];end
if nargin<11 || isempty(nX);nX=100;end
if nargin<12 || isempty(toler);toler=1e-8;end
if nargin<13 || isempty(lambda);lambda=0;end

%PRECONDITIONER
[P,SH]=precondDISORDER(S);

%MOTION FACTORS
if ~iscell(T)
    et=cell(1,2);
    if ~isempty(T)
        for l=1:2;et{l}=precomputeFactorsSincRigidTransform(kGrid,rkGrid,T,l-1,[],[],1);end
    end
else
    et=T;
end

if ~iscell(TFov)
    etFov=cell(1,2);
    if ~isempty(TFov)
        for l=1:2;etFov{l}=precomputeFactorsSincRigidTransform(kGrid,rkGrid,TFov,l-1,[],[],1);end
    end
else
    etFov=TFov;
end

n=0;
%DECODE THE RESIDUALS
if ~isempty(x)    
    r=decodeDISORDER(y,SH,et{1},etFov{1},A,B,M)-systemDISORDER(x);
    n=n+3;
else
    r=-decodeDISORDER(y,SH,et{1},etFov{1},[],B,M);
    n=n+1;
    NX=size(r);
    x=zeros(NX,'like',r);
end

z=P.*r;
p=z;
rsold=multDimSum(conj(z).*r,1:3);

%ITERATIONS
while 1
    Ap=systemDISORDER(p);
    n=n+2;
    al=conj(rsold)/multDimSum(conj(p).*Ap,1:3);
    xup=al*p;
    x=x+xup; 
    xup=max(abs(xup(:)).^2);
    if xup<toler || n>=nX;break;end
    r=r-al*Ap;
    z=P.*r;
    rsnew=multDimSum(conj(z).*r,1:3);
    be=rsnew/rsold;
    p=z+be*p;
    rsold=rsnew;
    if sqrt(abs(rsnew))<1e-10;break;end
    n=n+1;
end

if nargout>=3
    E=errorFit(x,y,S,et{2},A,B,kGrid,rkGrid);
    n=n+1;
end

function x=systemDISORDER(x)
    xS=encodeDISORDER(x,S,et{2},etFov{2},NY(1:2),B);
    xS=decodeDISORDER(xS,SH,et{1},etFov{1},A,B,M);
    x=xS+lambda*x;
end

end
