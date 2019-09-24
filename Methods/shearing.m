function x=shearing(x,k,dims)

%SHEARING shears the array x by k units along dimensions dims(1)-dims(2)
%   X=SHEARING(X,K,DIMS) translates the array X along dimension dims(1) 
%   with translations depending on the distance to the center on dimension 
%   dim(2)
%   * X is the array to be sheared
%   * K is the amount of shearing
%   * DIMS are the dimensions for shearing
%   ** X is the sheared array
%

gpu=isa(x,'gpuArray');
if any(k(:)~=0)
    N=size(x);N(setdiff(1:length(N),dims(2)))=1;
    rGrid=generateGrid(N,gpu,N,floor(N/2)+1);  
    H{dims(1)}=-bsxfun(@times,k,rGrid{dims(2)});
    x=shifting(x,H);
end
