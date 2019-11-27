function x=multDimMea(x,dim)

%MULTDIMMEA   Takes the mean of the elements of a multidimensional array 
%along a set of dimensions
%   X=MULTDIMMEA(X,{DIM})
%   * X is an array
%   * {DIM} are the dimensions over which to take the mean of the elements 
%   of the array, defaults to all
%   ** X is the contracted array
%

if nargin<2 || isempty(dim);dim=1:numDims(x);end

for n=1:length(dim);x=mean(x,dim(n));end
