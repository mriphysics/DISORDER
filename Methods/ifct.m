function x=ifct(x,n,dim)

%IFCT   Performs an Inverse Discrete Cosine Transform along a given 
%dimension. It is based on Andriy Myronenko's implementation in 
%https://uk.mathworks.com/matlabcentral/fileexchange/24050-multidimensional-discrete-cosine-transform--dct-
%   X=IFCT(X,{N},{DIM})
%   * X is the array to be transformed
%   * {N} is the transform size, so that X is padded with zeros if it has 
%   less than N points and truncated if it has more
%   * {DIM} is the dimension across which to take the transform
%   ** X is the transformed array
%

M=size(x);
if nargin<3 || isempty(dim);dim=find(M~=1,1);end
if isempty(dim);return;end

N=M;
if nargin>=2 && ~isempty(n);N(dim)=n;end

if N(dim)<M(dim);x=dynInd(x,1:N(dim),dim);elseif N(dim)>M(dim);x=padArray(x,N-M,0,'post');end
n=N(dim);

%PRECOMPUTE WEIGHTS   
if isa(x,'gpuArray')
    if ~isaUnderlying(x,'double');ww=single(2*exp((-1i*pi/(2*n))*(0:n-1)')/sqrt(2*n));else ww=2*exp((-1i*pi/(2*n))*(0:n-1)')/sqrt(2*n);end
else
    if ~isa(x,'double');ww=single(2*exp((-1i*pi/(2*n))*(0:n-1)')/sqrt(2*n));else ww=2*exp((-1i*pi/(2*n))*(0:n-1)')/sqrt(2*n);end
end
ww(1)=ww(1)/sqrt(2);
if isa(x,'gpuArray');ww=gpuArray(ww);end
perm=1:ndims(x);perm(1)=dim;perm(dim)=1;
ww=permute(ww,perm);

%INVERSE COSINE TRANSFORM
ind([1:2:n 2:2:n])=[1:ceil(n/2) n:-1:ceil(n/2)+1];
isrealx=isreal(x);
if ~isrealx;x=cat(length(M)+1,real(x),imag(x));end   
x=bsxfun(@times,ww,x);%Weight
x=fft(x,[],dim);%fft
x=real(dynInd(x,ind,dim));%Reorder
if ~isrealx;x=dynInd(x,1,length(M)+1)+1i*dynInd(x,2,length(M)+1);end

