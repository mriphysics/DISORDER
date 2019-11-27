function [z,w]=extractOrthogonalPlanes(x,y)

%EXTRACTORTHOGONALPLANES   Gets the orthogonal planes from a volume for 
%segmentation visualization
%   EXTRACTORTHOGONALPLANES(X,{Y})
%   * X is the image
%   * {Y} is the segmentation
%   * {PAU} indicates whether to pause the execution, it defaults to 1
%   * {ORTHV} shows orthogonal views
%   * {SL} serves to select a given slice
%   * {DY} serves to select a given dynamic
%

x=x(:,:,:,1);

N=size(x);
x=resampling(x,max(N)*ones(1,3),2);
N=size(x);
if ~isempty(y)
    y=y(:,:,:,1);
    y=resampling(y,max(N)*ones(1,3),2);    
    par=ellipsoidFromImage(single(y>0.5));
else
    par=mod(ceil(N/2)+1,N)+1;
end
z=zeros([N(1) N(2) 3],'like',x);
if ~isempty(y);w=zeros([N(1) N(2) 3],'like',y);end
for n=1:3
    x=shiftdim(x,1);
    z(:,:,n)=squeeze(dynInd(x,round(par(n)),2));
    if ~isempty(y)
        y=shiftdim(y,1);
        w(:,:,n)=squeeze(dynInd(y,round(par(n)),2));
    end
end
z=z(:,:);
if ~isempty(y);w=w(:,:);else w=[];end
