function x=multDimMin(x,dim)

%MULTDIMMIN   Takes the minimum of the elements of a multidimensional 
%array along a set of dimensions
%   X=MULTDIMMIN(X,{DIM})
%   * X is an array
%   * {DIM} are the dimensions over which to take the min of the elements 
%   of the array, defaults to all
%   ** X is the contracted array
%

if nargin<2 || isempty(dim);dim=1:numDims(x);end

for n=1:length(dim);x=min(x,[],dim(n));end