function [F,FH]=buildStandardDFTM(N,fn,gpu,re)

%BUILDSTANDARDDFTM   Builds standard DFT matrices
%   [F,FH]=BUILDSTANDARDDFTM(N,{FN},{GPU})
%   * N are the dimensions of the space
%   * {FN} indicates whether to generate fully unitary Fourier matrices. It
%   defaults to 0
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) matrices (empty, default depending on machine)
%   * {RE} is a flag to indicate that they are to be applied to a real 
%   image (defaults to 0)
%   ** F is a cell of discrete Fourier transform matrices along the 
%   different dimensions
%   ** FH is a cell of inverse discrete Fourier transform matrices along 
%   the different dimensions
%

if nargin<2 || isempty(fn);fn=0;end
if nargin<3 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
if nargin<4 || isempty(re);re=0;end

ND=length(N);
F=cell(1,ND);FH=cell(1,ND);
for m=1:ND
    if nargout>1;[F{m},FH{m}]=build1DFTM(N(m),fn,gpu,re);
    else F{m}=build1DFTM(N(m),fn,gpu,re);
    end
end
