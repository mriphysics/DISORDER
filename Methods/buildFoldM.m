function [A,AH]=buildFoldM(NS,NY,gpu,emre)

%BUILDFOLDM   Builds the fold matrices 
%   [A,AH]=BUILDFOLDM(NS,NY,{GPU},{EMRE})
%   * NS is the spatial size of the array
%   * NY is the spectral size of the array
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) matrices (empty, default depending on machine)
%   * {EMRE} serves to generate an empty return (it this feature is not 
%   to be used)
%   ** A is a cell of matrices to perform folding
%   ** AH is a cell of matrices to perform inverse folding
%

if nargin<3 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
if nargin<4 || isempty(emre);emre=0;end

ND=length(NS);
assert(ND==length(NY),'The dimensionality of the spatial space (%d) is not the same as for the spectral space (%d)',ND,length(NY));

A=cell(1,ND);AH=cell(1,ND);
if ~emre
    for n=1:ND
        oddFactSENSE=2*ceil((ceil(NS(n)/NY(n))-1)/2)+1;
        oFRed=(oddFactSENSE-3)/2;    
        oFRedT=oFRed*NY(n);
        over=(NS(n)-NY(n))/2;
        FOV=[floor(over) ceil(over)];
        if NS(n)==NY(n)
            A{n}=single(eye(NY(n)));
        else
            A{n}=single(zeros(NY(n),NS(n)));    
            A{n}(repmat(1:NY(n),[1 2*oFRed+1])+(FOV(2)-oFRed*NY(n):NS(n)-FOV(1)+oFRed*NY(n)-1)*NY(n))=1;
            A{n}([1:FOV(1)-oFRedT NY(n)-FOV(2)+1+oFRedT:NY(n)]+([NS(n)-FOV(1)+1+oFRedT:NS(n) 1:FOV(2)-oFRedT]-1)*NY(n))=1;
        end
        if gpu;A{n}=gpuArray(A{n});end    
        AH{n}=A{n}';
    end
end