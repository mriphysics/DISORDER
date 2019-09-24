function x=upsampling(x,dims,ups,odd)

%UPSAMPLING upsamples the array x
%   X=UPSAMPLING(X,{DIMS},{UPS},{ODD}) upsamples the array X by
%   introducing UPS zeros in the even samples along dimensions DIMS
%   * X is the array to be upsampled
%   * {DIMS} are the dimensions for upsampling. It defaults to all
%   * {UPS} is the upsampling factor. It defaults to 2
%   * {ODD} indicates whether the resulting array is odd or even. It 
%   defaults to 1 (odd).
%   ** X is the upsampled array
%

ND=numDims(x);
if nargin<2 || isempty(dims);dims=1:ND;end
if nargin<3 || isempty(ups);ups=2;end
if nargin<4 || isempty(odd);odd=1;end
N=size(x);N(end+1:max(dims))=1;

for n=1:length(dims)
    N(dims(n))=ups*N(dims(n))-(ups-1)*odd;
    xups=zeros(N,'like',x);
    xups=dynInd(xups,1:ups:N(dims(n)),dims(n),x);
    x=xups;
end
