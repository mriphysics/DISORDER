function x=mapMat(x,dims,col)

%MAPMAT   Maps a multidimensional array onto a matrix 
%   X=MAPMAT(X,DIMS,{COL})
%   * X is an array
%   * DIMS are a set of dimensions to be mapped onto the rows/columns; they
%   will be mapped in the order specified by dims
%   * {COL} indicates whether the dimensions have to be mapped to the
%   columns (1, default) or to the rows (0)
%   ** X is the reshaped result with fully preserved order of dimensions
%   other than the mapped
%

if nargin<3 || isempty(col);col=1;end

if isempty(x);return;end
assert(~isempty(dims) && dims>=1,'Dimensions to be mapped have to be meaningful');

Nper=length(dims);
N=size(x);N(end+1:end+Nper)=1;

nodims=1:ndims(x)+Nper;nodims(dims)=[];
if ~col
    if max(dims)~=1
        perm(1:Nper)=dims;perm(Nper+1:ndims(x)+Nper)=nodims;
        x=permute(x,perm);
    end
    x=reshape(x,[prod(N(dims)) prod(N(nodims))]);
else
    perm(1:ndims(x))=nodims;perm(ndims(x)+1:ndims(x)+Nper)=dims;
    x=permute(x,perm);
    x=reshape(x,[prod(N(nodims)) prod(N(dims))]);
end