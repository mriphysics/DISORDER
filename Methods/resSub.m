function [x,N]=resSub(x,orDim,M)

%RESSUB   Reshapes a set of dimensions of an array and extends/reduces the 
%array accordingly
%   [X,N]=RESSUB(X,ORDIM,{M})
%   * X is an array
%   * ORDIM are a set of contiguous dimensions that need to be rearranged
%   * {M} is the size of the new set of dimensions to be created. It should
%   multiply to the size of dimensions given by ORDIM. Defaults to a
%   singleton with the product of the sizes of dimensions given by ORDIM
%   ** X is the reshaped and populated result
%   ** N is the new size of X
%

Ndo=length(orDim);
if Ndo>0
    %assert(all(orDim>0),'Origin dimensions have to be positive');%THIS IS COMMENTED FOR ACCELERATED RUNS

    N=size(x);N(end+1:max(orDim))=1;
    if nargin<3 || isempty(M);M=prod(N(orDim));end

    %assert(prod(M)==prod(N(orDim)),'Size of reshaping (%d) has to match size of original dimensions (%d)',prod(M),prod(N(orDim)));% THIS IS COMMENTED FOR ACCELERATED RUNS
    %assert(~(Ndo>1 && any(diff(orDim)~=1)),'The origin dimensions (%s) have to be contiguous',sprintf('%d ',orDim));%THIS IS COMMENTED FOR ACCELERATED RUNS

    dimN=[M N(orDim(end)+1:end)];
    N(orDim(1):orDim(1)+length(dimN)-1)=dimN;
    N(orDim(1)+length(dimN):end)=1;
    x=reshape(x,N);
end
N=size(x);