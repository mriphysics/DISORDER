function x=plugNoise(x,useR)

% PLUGNOISE plugs complex circularly symmetric AWGN noise samples in an 
%array
%   X=PLUGNOISE(X,USER) 
%   * X is the input array
%   * USER enables the usage of real only noise
%   ** X is the noise array with same dimensions as the input array
%

gpu=isa(x,'gpuArray');
if nargin<2;useR=0;end%Complex only, for back-compatibility!

comp=~isreal(x) || ~useR;

N=size(x);
if comp;x=single(randn(N))+1i*single(randn(N));else x=single(randn(N));end
if gpu;x=gpuArray(x);end
