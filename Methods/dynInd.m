function x=dynInd(x,ind,dim,y)

%DYNIND   Performs dynamic indexing over multidimensional arrays
%   X=DYNIND(X,IND,DIM,{Y})
%   * X is an array
%   * IND is a vector or a cell of vectors of length equal to the length of
%   DIM. If a vector and DIM a singleton it extracts the indexes in the 
%   vector from array X. If a vector and DIM non singleton it has to have 
%   the same dimensions as DIM and will extract the corresponding index 
%   from the dimensions given by DIM. If a cell it extracts the indexes in 
%   each element of the cell from the corresponding dimension in DIM
%   * DIM are the dimensions for indexing
%   * {Y} is a second input from which all its elements are assigned to the
%   corresponding indexes in X
%   ** X is the indexed result
%
%   A bit about sintaxis:
%   We create a multidimensional array
%   x=randn(5,5,5,5,5);
%   Vector mode of operation along a single dimension
%   size(dynInd(x,[3 2],3))
%
%    ans =
%
%         5     5     2     5     5
%
%   Vector mode of operation along several dimensions (only scalar indexing)
%   size(dynInd(x,[3 2],[3 4]))
%
%    ans =
%
%         5     5     1     1     5
%
%   Cell mode of operation along a single dimension, equivalent to the 
%   vector mode of operation along a single dimension
%   size(dynInd(x,{[3 2]},3))
%
%    ans =
%
%         5     5     2     5     5
%
%   Cell mode of operation along several dimensions, more flexible than 
%   vector mode of operation along several dimensions as it allows 
%   vectorial indexing
%   size(dynInd(x,{[3 2],1},[3 4]))
%
%   ans =
%
%        5     5     2     1     5
%

Ndim=length(dim);
%if iscell(ind);assert(length(ind)==Ndim,'Lenght of cell of indexes (%d) has to match length of dimensions to index (%d)',length(ind),Ndim);else assert(~(Ndim>1 && length(ind)~=Ndim),'Length of vector of indexes (%d) has to match the length of dimensions to index (%d). You may want to consider operating with cells of indexes for vectorial indexing',length(ind),Ndim);end%Takes time and may simply work without error control
%assert(length(unique(dim))==length(dim),'Dimensions to index do not allow repetitions');%Takes time and may simply work without error control

ndx=max(ndims(x),max(dim));

if ndx==1;subs={':'};%Calling repmat takes time
elseif ndx==2;subs=[{':'} {':'}];
elseif ndx==3;subs=[{':'} {':'} {':'}];
elseif ndx==4;subs=[{':'} {':'} {':'} {':'}];
elseif ndx==5;subs=[{':'} {':'} {':'} {':'} {':'}];
elseif ndx==6;subs=[{':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==7;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==8;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==9;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==10;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==11;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==12;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==13;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==14;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
else subs=repmat({':'}, [1 ndx]);
end

%Note error control is not implemented, if the user tries to index outside 
%the size of x along a given dimension, matlab should throw an error
if ~iscell(ind)
    if Ndim~=1;subs(dim)=num2cell(ind);else subs{dim}=ind;end
else
    subs(dim)=ind;%If this fails it may be because the length of the cell does not match the length of dimensions
end
if nargin<4;x=x(subs{:});else x(subs{:})=y;end

