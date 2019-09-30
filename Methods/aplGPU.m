function x=aplGPU(A,x,m)

%APLGPU   Applies the matrix A to the array x along dimension m of x
%   X=APLGPU(A,X,M)
%   * A is a matrix to apply over A.
%   * X is the array on which to apply A
%   * M is the direction along which to apply A over X
%   ** X is the transformed array
%

if isempty(x) || isempty(A);return;end

assert(size(A,2)==size(x,m),'Number of columns of the matrix (%d) does not match array size (%d) for dim %d',size(A,2),size(x,m),m);
[M,N]=parUnaFun({A,x},@size);
assert(ismatrix(A),'Operator should be a matrix');
N(max(end+1,m+1):4)=1;

if isscalar(A)
    x=A*x;
else
    if m==1;x=x(:,:);else x=mapMat(x,m,0);end
    x=A*x;
    if m==1;x=reshape(x,[M(1) N(2:end)]);else x=resPop(x,[1 2],[M(1) N(1:m-1) N(m+1:length(N))],[m 1:m-1 m+1:length(N)]);end
end

