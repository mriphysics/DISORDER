function x=multDimMax(x,dim)

%MULTDIMMAX   Takes the maximum of the elements of a multidimensional 
%array along a set of dimensions
%   X=MULTDIMMAX(X,{DIM})
%   * X is an array
%   * {DIM} are the dimensions over which to take the max of the elements of 
%   the array, defaults to all
%   ** X is the contracted array
%

if nargin<2 || isempty(dim);dim=1:numDims(x);end

for n=1:length(dim);x=max(x,[],dim(n));end