function [rGrid,kGrid,rkGrid,kkGrid,cGrid]=generateTransformGrids(N,gpu,Nres,cent,sh,prods,NSlab,slabInd)

%GENERATETRANSFORMGRIDS   Generates spatial/k-spatial grids to apply a 
%rigid transform
%   [RGRID,KGRID,RKGRID]=GENERATETRANSFORMGRIDS(N,{GPU},{NRES},{CENT},{SH},{PRODS},{NSLAB},{SLABIND})
%   * N are the dimensions of the space
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) arrays (empty, default depending on machine)
%   * {NRES} defines the size of the output grid (defaults to N)
%   * {CENT} defines the center of coordinates (defaults to the center of 
%   the FOV)
%   * {SH} is a flag to indicate whether to shift the grids (defaults to 0)
%   * {PRODS} is an array of values to build grids with multiple sizes
%   * {NSLAB} serves to use slab-based transforms
%   * {SLABIND} extracts a series of slices around which to build the slabs
%   ** RGRID is the spatial grid (from -N/2 to N/2)
%   ** KGRID is the spectral grid (from -pi to pi)
%   ** RKGRID is the spatio-spectral grid
%   ** KKGRID is the spectral-spectral grid
%   ** CGRID is the center of the grid
%

N(end+1:3)=1;N=N(1:3);
if nargin<2 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
if nargin<3 || isempty(Nres);Nres=round(N);end
if nargin<4 || isempty(cent);cent=(Nres/2)+1;end
if nargin<5 || isempty(sh);sh=0;end
if nargin<6;prods=[];end
if nargin<7;NSlab=[];end
if nargin<8;slabInd=[];end

if ~isempty(NSlab)    
    rGrid=generateTransformGrids(N,gpu,Nres,cent,sh,prods);
    NResAux=N;NResAux(3)=NSlab;NAux=N;NAux(3)=NSlab;    
    [~,kGrid]=generateTransformGrids(NAux,gpu,NResAux,[],sh);
    rGrid{3}=extractSlabs(rGrid{3},NSlab,1,1);    
    if ~isempty(slabInd);rGrid{3}=dynInd(rGrid{3},slabInd,6);end
    cGrid=[];
else
    if ~isempty(prods)
        NNres=N./Nres;
        if all(ceil((Nres+1)/2)==cent);flcentNres=1;else flcentNres=0;centN=cent./Nres;end
        for n=1:3;Nres(n)=nextprod(prods,Nres(n));end    
        N=Nres.*NNres;
        if ~flcentNres;cent=Nres.*centN;else cent=ceil((Nres+1)/2);end
    end

    [rGrid,~,~,cGrid]=generateGrid(Nres,gpu,N,cent);
    kGrid=generateGrid(Nres,gpu,2*pi*Nres./N,ceil((Nres+1)/2));
    if sh
        for m=1:3;kGrid{m}=ifftshift(kGrid{m},m);end
    end
end

per=[1 3 2;
     2 1 3];
rkGrid=cell(2,3);
for n=1:2
    for m=1:3
        rkGrid{n}{m}=bsxfun(@times,rGrid{per(3-n,m)},kGrid{per(n,m)});
    end
end

fact=[1 2 3 1 1 2;
      1 2 3 2 3 3];
kkGrid=cell(1,6);
for m=1:6;kkGrid{m}=bsxfun(@times,kGrid{fact(1,m)},kGrid{fact(2,m)});end