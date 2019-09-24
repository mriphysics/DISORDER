function nd=numDims(x)

%NUMDIMS   Modifies the default behaviour of the ndims function for empty 
%arrays and row-like arrays, allowing to return a number of dimensions 
%lower than 2
%   ND=NUMDIMS(X)
%   * X is an array
%   ** ND are the dimensions of the array
%

if isempty(x);nd=-1;return;
elseif numel(x)==1;nd=0;return;
elseif size(x,1)>1 && size(x,2)==1 && ismatrix(x);nd=1;return;
else nd=ndims(x);
end