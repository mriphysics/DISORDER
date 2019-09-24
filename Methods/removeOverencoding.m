function x=removeOverencoding(x,overEnc)

% REMOVEOVERENCODING removes excess FOV
%   X=REMOVEOVERENCODING(X,OVERDEC)
%   * X is the data to be cropped
%   * OVERENC is the overencoding factor, see variable Alg.OverDec in file 
%   reconAlgorithm.m
%   ** X is the cropped data
%

overEnc(overEnc<0)=1;%To not remove if negative in which case we write all the information
N=size(x);N(end+1:3)=1;N=N(1:length(overEnc));
x=resampling(x,round(N./overEnc),2);
