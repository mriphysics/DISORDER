function writeNII(dat,suffix,x,MS,MT)

%WRITENII   Writes a set of nii files
%   WRITENII(DAT,SUFF,X,{MS},{MT})
%   * DAT indicates the path and name of the file
%   * SUFF is a list of suffixes of the files to be written
%   * X is the data to write
%   * {MS} is the spacing of the data, defaults to 1
%   * {MT} is the orientation of the data, defaults to identity
%

%ERROR CONTROL AND DEFAULT STRUCTURES
N=length(suffix);
if ~exist('MS','var')
    MSaux=[1 1 1];
    MS=cell(1,N);for n=1:N;MS{n}=MSaux;end
end
if ~exist('MT','var')
    MTaux=eye(4,4);    
    MT=cell(1,N);
    for n=1:N
        for m=1:3;MTaux(m,m)=MS{n}(m);end
        MT{n}=MTaux;
    end
end
NIn=[length(suffix) length(x) length(MS) length(MT)];
assert(length(unique(NIn))==1,'Size of input cells must be the same. Currently it is: %d %d %d %d',length(suffix),length(x),length(MS),length(MT));

%WRITE IMAGES
for n=1:N    
    M=size(x{n});M(end+1:12)=1;
    x{n}=reshape(single(x{n}),[M(1:3) prod(M(4:12))]);    
    if any(imag(x{n}(:)))
        p=gather(angle(x{n}));       
        niftiIm=make_nii(p,MS{n}(1:3));
        niftiIm.hdr.hist.srow_x=MT{n}(1,:);niftiIm.hdr.hist.srow_y=MT{n}(2,:);niftiIm.hdr.hist.srow_z=MT{n}(3,:);niftiIm.hdr.hist.sform_code=1;         
        if length(MS{n})>=4;niftiIm.hdr.dime.pixdim(5)=MS{n}(4);end
        save_nii(niftiIm,sprintf('%s_%sPh.nii',dat,suffix{n}));
        x{n}=abs(x{n});
    end
    x{n}=gather(real(x{n}));    
    niftiIm=make_nii(x{n},MS{n}(1:3));
    niftiIm.hdr.hist.srow_x=MT{n}(1,:);niftiIm.hdr.hist.srow_y=MT{n}(2,:);niftiIm.hdr.hist.srow_z=MT{n}(3,:);niftiIm.hdr.hist.sform_code=1;            
    if length(MS{n})>=4;niftiIm.hdr.dime.pixdim(5)=MS{n}(4);end
    save_nii(niftiIm,sprintf('%s_%s.nii',dat,suffix{n}));
end
