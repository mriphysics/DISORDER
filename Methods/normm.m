function x=normm(x,y,dim)

%NORMM   Computes the squared l2-norm of an array or the squared 
%l2-distance between a pair of arrays
%   X=NORMM(X,{Y},{DIM})
%   * X is the input array
%   * {Y} is another array (to compute differences)
%   * {DIM} serves to compute summation only along certain dimensions
%   ** X is the squared norm
%

if nargin>=2 && ~isempty(y);x=bsxfun(@minus,x,y);end
ND=numDims(x);
if nargin<3 || isempty(dim);dim=1:ND;end
x=multDimSum(abs(x).^2,dim);

