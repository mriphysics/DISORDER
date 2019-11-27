function [pGrid,DeltaX,DeltaF,cent]=generateGrid(N,gpu,nor,cent,shift)

%GENERATEGRID   By default it generates a normalized grid of size N
%spanning the range [-0.5 0.5]
%   [PGRID,DELTAX,DELTAF,CENT]=GENERATEGRID(N,{GPU},{NOR},{CENT},{SHIFT})
%   * N are the dimensions of the space
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) arrays (empty, default depending on machine)
%   * {NOR} normalizes the grid so that it spans the space [-NOR/2 NOR/2]. It
%   defaults to 1. Normalization can be different for different dimensions
%   * {CENT} centers the grid at the corresponding discrete point. Defaults
%   to (N+1)/2. Note that for even N the center may not be included in the
%   grid. To do so, one could opt for a ceil((N+1)/2) for DFT grid types
%   * {SHIFT} on top of previous basic parameters, it allows to shift the
%   location of the grid by a given number of discrete points.
%   ** PGRID is the generated grid
%   ** DELTAX is the resolution of the grid
%   ** DELTAF is the FOV of the grid
%   ** CENT is the center of the grid
%

if nargin<2 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
if nargin<3 || isempty(nor);nor=1;end
if nargin<4 || isempty(cent);cent=(N+1)/2;end
if nargin<5 || isempty(shift);shift=0;end

cent=cent-shift;
NDims=length(N);
if length(cent)~=NDims;fprintf('The center of the grid has size %d while the grid has size %d\n',length(cent),NDims);end
pGrid=cell(1,NDims);
if length(nor)==1;nor=repmat(nor,[1 NDims]);end
DeltaX=single(ones(1,NDims));
for m=1:NDims
    pGrid{m}=single((nor(m)/N(m))*((1:N(m))-cent(m)));
    if gpu>0;pGrid{m}=gpuArray(pGrid{m});end
    DeltaX(m)=nor(m)/N(m);
    perm=1:NDims;perm(2)=m;perm(m)=2;
    pGrid{m}=permute(pGrid{m},perm);
end
DeltaF=DeltaX.*N;
