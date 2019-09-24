function [x,MS,MT]=readNII(dat,suff,gpu)

%READNII   Reads a set of nii files
%   [X,MS,MT]=READNII(DAT,SUFF,{GPU})
%   * DAT indicates the path and name of the file
%   * SUFF is a list of suffixes of the files to be read
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu 
%   (0) arrays (empty, default depending on machine)
%   * X is the returned data
%   * MS is the returned spacing
%   * MT is the returned orientation
%

if ~exist('gpu','var') || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end

N=length(suff);
x=cell(1,N);MS=cell(1,N);MT=cell(1,N);
for n=1:N
    nii=load_untouch_nii(sprintf('%s_%s.nii',dat,suff{n}));
    x{n}=nii.img;
    if gpu;x{n}=gpuArray(x{n});end
    if exist(sprintf('%s_%sPh.nii',dat,suff{n}),'file')
        nii=load_untouch_nii(sprintf('%s_%sPh.nii',dat,suff{n}));
        p=nii.img;
        if gpu;p=gpuArray(p);end
        x{n}=x{n}.*exp(1i*p);
    end
    MS{n}=nii.hdr.dime.pixdim(2:4);
    MT{n}=eye(4);MT{n}(1,:)=nii.hdr.hist.srow_x;MT{n}(2,:)=nii.hdr.hist.srow_y;MT{n}(3,:)=nii.hdr.hist.srow_z;
end
