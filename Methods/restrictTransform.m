function T=restrictTransform(T,limTran,ndT)

%RESTRICTTRANSFORM   Restricts the angular components of the transform to 
%   lie in the range -pi/pi
%   T=RESTRICTTRANSFORM(T,{LIMTRAN},{NDT})
%   * T is the transform to be restricted
%   * {LIMTRAN} are limits for translation, it defaults to emtpy
%   * {NDT} is the dimension with the transform parameters
%   ** T is the restricted transform
%

if nargin<3 || isempty(ndT);ndT=ndims(T);end

Tang=dynInd(T,4:6,ndT);
if nargin>=2
    perm=1:ndT;perm([1 ndT])=[ndT 1];
    limAng=[-pi pi;-pi pi;-pi pi];
    limAng=permute(limAng,perm);  
    Tang(bsxfun(@ge,Tang,dynInd(limAng,2,2)))=0;
    Tang(bsxfun(@le,Tang,dynInd(limAng,1,2)))=0;
    
    limTran=permute(limTran,perm);
    Ttra=dynInd(T,1:3,ndT);
    Ttra(bsxfun(@ge,Ttra,dynInd(limTran,2,2)))=0;
    Ttra(bsxfun(@le,Ttra,dynInd(limTran,1,2)))=0;
    T=dynInd(T,1:3,ndT,Ttra);
else
    while ~isempty(Tang(Tang>pi));Tang(Tang>pi)=Tang(Tang>pi)-2*pi;end
    while ~isempty(Tang(Tang<-pi));Tang(Tang<-pi)=Tang(Tang<-pi)+2*pi;end
end
    
T=dynInd(T,4:6,ndT,Tang);