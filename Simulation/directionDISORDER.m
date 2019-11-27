function [Tup,w,flagw,winic]=directionDISORDER(dt,ddt,w,flagw,winic)

%DIRECTIONDISORDER   Computes the motion update for the motion parameters
%   [TUP,W,DTP,STP]=DIRECTIONDISORDER(DT,DDT,W,DTP,STP)
%   * DT is the gradient direction
%   * DDT is the Hessian
%   * W is the update weight
%   * FLAGW flags current energy decrease status
%   * WINIC is the baseline LM damping parameter
%   ** TUP is the motion update
%   ** W is the update weight
%   ** WINIC is the baseline LM damping parameter
%

maxTDeg=1*pi/180;%Maximum motion update in degrees (for CG-like)
multA=1.2;%Factor to divide the weight when E(end)<=E(end-1)
multB=2;%Factor to multiplicate the weight when H(end)>H(end-1)

ND=size(dt);
NDD=size(ddt);

%WEIGHT INITIALIZATION
if isempty(w)
    if isempty(winic);w=repmat(1e-3,[1 ND(2)]);else w=repmat(winic,[1 ND(2)]);end
end%LM-like

if isempty(flagw)
    flagw=zeros([1 ND(2)]);
else
    %WEIGHT UPDATE
    w(flagw==1)=w(flagw==1)/multA;%Energy decrease
    w(flagw==0)=w(flagw==0)*multB;%No energy decrease    
end

ha=[1 2 3 1 1 2 1 2 3 1 2 3 1 2 3 4 4 5 4 5 6;
    1 2 3 2 3 3 4 4 4 5 5 5 6 6 6 5 6 6 4 5 6];%Second order derivatives ordered as in sincRigidTransformHessian

ddtd=repmat(eye(ND(1)),[1 1 ND(2)]);%Damped Hessian
w=permute(w,[1 3 2]);
dt=permute(dt,[1 3 2]);
ddt=permute(ddt,[1 3 2]);
for k=1:NDD(1)
    if ha(1,k)==ha(2,k)
        ddtd(ha(1,k),ha(2,k),:)=bsxfun(@times,ddt(k,1,:),(1+w));                       
    else
       ddtd(ha(1,k),ha(2,k),:)=ddt(k,1,:);
       ddtd(ha(2,k),ha(1,k),:)=ddt(k,1,:);
    end
end 
ddtd=ddtd+bsxfun(@times,multDimMax(abs(ddtd),1:2)/1e6,eye(ND(1)));%To stabilize
st=single(matfun(@mldivide,double(ddtd),double(dt)));

dH=bsxfun(@rdivide,st,w);
if isempty(winic);winic=maxTDeg/multDimMax(abs(dH(4:6,:,:)),1:3);end    
dH=winic*dH;

Tup=permute(dH,[2 4 5 6 3 1]);
flagw(:)=0;
w=permute(w,[1 3 2]);
