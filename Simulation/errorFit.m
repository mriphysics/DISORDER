function E=errorFit(x,y,S,T,TFov,A,B,kGrid,rkGrid,sep,M)

%ERRORFIT   Measures the error for the least-squares reconstruction
%formulation
%   E=ERRORFIT(X,Y,S,{T},{TFOV},{A},{B},{KGRID},{RKGRID},{SEP},{M}) returns the value of the 
%   objective function of the reconstruction
%   * X is the reconstructed image
%   * Y is the measured data
%   * S is the coil-array sensitivity map
%   * {T} are the estimated transforms for each shot
%   * {TFOV} are the transforms with the geometry for each repetition
%   * {A} is the sampling scheme
%   * {B} is the multiband sampling scheme
%   * {KGRID} is a grid of points in the spectral domain
%   * {RKGRID} is a grid of points in the spatial-spectral domain
%   * {SEP} does not perform summation along the specified dimensions
%   * {M} serves to extract certain harmonic components
%   ** E is the value of the residuals
%

NY=size(y);
if nargin<4;T=[];end
if nargin<5;TFov=[];end
if nargin<6;A=NY(1:2);end
if nargin<7 || isempty(B);B=cell(1,2);end
if nargin<8;kGrid=[];end
if nargin<9;rkGrid=[];end
if nargin<10;sep=[];end
if nargin<11;M=[];end

if ~iscell(T)
    if ~isempty(T);et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,T,1,[],[],1);else et=[];end
else
    et=T;
end
if ~iscell(TFov)
    if ~isempty(TFov);etFov=precomputeFactorsSincRigidTransform(kGrid,rkGrid,TFov,1,[],[],1);else etFov=[];end
else
    etFov=TFov;
end
x=encodeDISORDER(x,S,et,etFov,A,B,y);
if ~isempty(M);x=bsxfun(@times,M,x);end
E=normm(x,[],setdiff(1:8,sep));

