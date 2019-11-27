function [dt,ddt]=gradientDISORDER(Jt,ry,S,etFov,A,B)

%GRADIENTDISORDER   Computes the gradient of with respect to the motion
%parameters
%   [DT,DDT]=GRADIENTDISORDER(JT,RY,{S},{ETFOV},{A},{B})  
%   * JT is the Jacobian with respect to the motion parameters
%   * RY are the residuals
%   * {S} is the coil array sensitivity map
%   * {ETFOV} are the parameters of the transform to the acquired geometry 
%   space
%   * {A} contains the sampling scheme or the size of the Fourier space
%   * {B} contains the multiband sampling scheme
%   ** DT is the returned gradient
%   ** DDT is the returned Hessian
%

NY=size(ry);

if nargin<3;S=[];end
if nargin<4;etFov=[];end
if nargin<5;A=NY(1:2);end
if nargin<6 || isempty(B);B=cell(1,2);end

ha=[1 2 3 1 1 2 1 2 3 1 2 3 1 2 3 4 4 5 4 5 6;
    1 2 3 2 3 3 4 4 4 5 5 5 6 6 6 5 6 6 4 5 6];%Second order derivatives ordered as in sincRigidTransformHessian
ND=[6 size(Jt{1},5)*size(Jt{1},6)*size(Jt{1},7)*size(Jt{1},8)];
NDD=[21 ND(2)];
dt=zeros(ND);ddt=zeros(NDD);
for m=1:ND(1)       
    if NY(3)~=1 || ismember(m,[1 2 4]);Jt{m}=encodeDISORDER(Jt{m},S,[],etFov,A,B);end
end
for m=1:NDD(1)
    if NY(3)~=1 || all(ismember(ha(:,m),[1 2 4]));ddt(m,:)=resSub(permute(gather(multDimSum(real(Jt{ha(1,m)}.*conj(Jt{ha(2,m)})),1:4)),[1 5 6 7 8 3 2 4]),2:5);end
end
for m=1:ND(1)
    if NY(3)~=1 || ismember(m,[1 2 4]);dt(m,:)=resSub(permute(gather(multDimSum(real(bsxfun(@times,Jt{m},conj(ry))),1:4)),[1 5 6 7 8 3 2 4]),2:5);end
end

