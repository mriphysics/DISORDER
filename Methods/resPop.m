function [x,N]=resPop(x,orDim,M,deDim)

%RESPOP   Reshapes a set of dimensions of an array and populates another 
%   set of dimensions with the result
%   [X,N]=RESPOP(X,ORDIM,{M},{DEDIM})
%   * X is an array
%   * ORDIM are a set of dimensions that need to be rearranged
%   * {M} is the size of the new set of dimensions to be created. It should
%   multiply to the size of dimensions given by ORDIM. If empty it 
%   corresponds to the size of the input array for dimension ORDIM
%   * {DEDIM} is the set of dimensions where to place the results. If
%   ommited or empty the dimensions are moved to the end of the array
%   ** X is the reshaped and populated result
%   ** N is the new size of X
%

N=size(x);NDor=ndims(x);

if nargin<3 || isempty(M)
    M=N(orDim(orDim<=NDor));M=[M ones(sum(orDim>NDor))];
    if length(deDim)==1;M=prod(M);end
end
if nargin<4;deDim=[];end

indTrim=find(M~=1,1,'last');
if ~isempty(indTrim) && indTrim~=1
    M(indTrim+1:end)=[];
    if ~isempty(deDim);deDim(indTrim+1:end)=[];end
end

[Ndr,Ndo]=parUnaFun({M,orDim},@length);
if isempty(deDim);deDim=NDor+1:NDor+Ndr;end
assert(length(deDim)==Ndr,'Number of destination dimensions (%d) has to match dimensionality of reshaping (%d)',length(deDim),Ndr);
assert(all(orDim>0),'Origin dimensions have to be positive');
assert(all(deDim>0),'Destination dimensions have to be positive');

N(end+1:max(max(end+1,max(orDim)),max(deDim)))=1;
assert(prod(M)==prod(N(orDim)),'Size of reshaping (%d) has to match size of original dimensions (%d)',prod(M),prod(N(orDim)));
assert(~(Ndo>1 && any(diff(orDim)~=1)),'The origin dimensions (%s) have to be contiguous',sprintf('%d ',orDim));

%assert(all(N(deDim(~ismember(deDim,orDim)))==1),'The destination dimensions (%s) have to be all singletons and not of size (%s)',sprintf('%d ',deDim),sprintf('%d ',N(deDim)));%This takes a bit of time, commented, hope it does not break...

ND=length(N);

perm=1:(ND+Ndr);perm(orDim)=Ndr+(ND+1:ND+Ndo);perm(Ndr+(ND+1:ND+Ndo))=orDim;
x=permute(x,perm);
L=[N(1:min(orDim)-1) ones(1,Ndo) N(max(orDim)+1:ND) M];
x=reshape(x,L);
NOL=length(L);
perm=1:NOL;perm(deDim)=NOL-Ndr+1:NOL;perm(NOL-Ndr+1:NOL)=deDim;
x=permute(x,perm);
N=size(x);
