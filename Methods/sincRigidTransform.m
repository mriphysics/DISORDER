function [x,xB]=sincRigidTransform(x,et,di,F,FH,sumR,useFull,re,padfl)

%SINCRIGIDTRANSFORM rigidly transforms volumes using sinc-based interpolation (both forwards and backwards)
%   [X,XB]=SINCRIGIDTRANSFORM(X,ET,{DI},{F},{FH},{SUMR},{USEFULL})
%   * X is a volume
%   * ET are the transform factors
%   * {DI} is a flag to indicate whether to perform direct (default) or 
%   inverse transform
%   * {F} contains discrete Fourier transform matrices
%   * {FH} contains inverse discrete Fourier transform matrices
%   * {SUMR} indicates whether to sum along the motion states when 
%   applying the inverse transform (defaults to 1)
%   * {USEFULL} indicates whether to use the full transform (1, default),
%   only the translations (0) or only the rotations (-1)
%   * {RE} is a flag to indicate that the transforms are to be applied to a
%   real image (defaults to 0). THIS IS DEPRECATED IT SHOULD ALWAYS BE 0
%   * {PAD} is a flag to indicate that the transforms have to be padded 
%   (defaults to 0)
%   ** X is the rigidly transformed volume
%   ** XB are intermediate results that can be reused to compute the gradient of the transform
%

gpu=isa(x,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

if nargin<3 || isempty(di);di=1;end
if nargin<4 || isempty(F);F=cell(1,3);end
if nargin<5 || isempty(FH);FH=cell(1,3);end
if nargin<6 || isempty(sumR);sumR=1;end
if nargin<7 || isempty(useFull);useFull=1;end
if nargin<8 || isempty(re);re=0;end
if nargin<9 || isempty(padfl);padfl=0;end

if iscell(et{1});for n=1:3;Npad(n)=size(et{1}{n},n);end
else Npad=size(et{1});Npad(end+1:3)=1;Npad=Npad(1:3);
end
N=size(x);N(end+1:3)=1;N=N(1:3);
x=resampling(x,Npad,2);

if padfl && gpuF==2;gpuF=1;end

tr=[1 3 2;
    2 1 3];

if ~iscell(et{1});ndT=max(ndims(et{1}),4);else ndT=max(ndims(et{1}{1}),4);end

if di==1
    xB=cell(1,5);xB{5}=x;        
    %ROTATIONS
    if useFull~=0
        for m=1:3  %Axis: 3-2-1: for LR(fast PE)-AP(slow PE)-FH(read) it would be yaw-roll-pitch
            if any(et{5}{m}(:)==1);x=bsxfun(@times,x,1-et{5}{m})+bsxfun(@times,flipping(x,et{4}{m}),et{5}{m});end%FLIPPING FOR LARGER THAN 90DEG ROTATIONS
            x=fftGPU(x,tr(1,m),gpuF,F{tr(1,m)});
            xB{m}=x;    
            x=bsxfun(@times,x,et{2}{m});
            x=ifftGPU(x,tr(1,m),gpuF,FH{tr(1,m)});
            x=fftGPU(x,tr(2,m),gpuF,F{tr(2,m)});
            x=bsxfun(@times,x,et{3}{m});
            x=ifftGPU(x,tr(2,m),gpuF,FH{tr(2,m)});
            x=fftGPU(x,tr(1,m),gpuF,F{tr(1,m)});
            x=bsxfun(@times,x,et{2}{m});
            if m~=3 || re;x=ifftGPU(x,tr(1,m),gpuF,FH{tr(1,m)});end
        end
    end
    %TRANSLATION
    if useFull~=-1
        if ~iscell(et{1})
            for m=1:3     
                if m~=tr(1,3) || useFull==0 || re;x=fftGPU(x,m,gpuF,F{m});end             
            end
            xB{4}=x;
            x=bsxfun(@times,x,et{1});     
            for m=3:-1:1;x=ifftGPU(x,m,gpuF,FH{m});end
        else
            for m=1:3
                if m~=tr(1,3) || useFull==0 || re;x=fftGPU(x,m,gpuF,F{m});end
            end
            xB{4}=x;       
            for m=1:3;x=bsxfun(@times,x,et{1}{m});end
            for m=3:-1:1;x=ifftGPU(x,m,gpuF,FH{m});end     
        end            
    elseif ~re
        x=ifftGPU(x,tr(1,3),gpuF,FH{tr(1,3)});
    end
else
    %BACK-TRANSLATION
    if useFull~=-1
        if ~iscell(et{1})
            for m=1:3;x=fftGPU(x,m,gpuF,F{m});end
            x=bsxfun(@times,x,et{1});         
            for m=3:-1:1  
                if m~=tr(1,3) || useFull==0 || re;x=ifftGPU(x,m,gpuF,FH{m});end
            end
        else
            for m=1:3
                x=fftGPU(x,m,gpuF,F{m});
                x=bsxfun(@times,x,et{1}{m});
                if m~=tr(1,3) || useFull==0 || re;x=ifftGPU(x,m,gpuF,FH{m});end
            end
        end
            
    elseif ~re
        x=fftGPU(x,tr(1,3),gpuF,F{tr(1,3)});
    end
    
    %BACK-ROTATION
    if useFull~=0
        for m=3:-1:1 
            if m~=3 || re;x=fftGPU(x,tr(1,m),gpuF,F{tr(1,m)});end
            x=bsxfun(@times,x,et{2}{m});
            x=ifftGPU(x,tr(1,m),gpuF,FH{tr(1,m)});
            x=fftGPU(x,tr(2,m),gpuF,F{tr(2,m)});
            x=bsxfun(@times,x,et{3}{m});
            x=ifftGPU(x,tr(2,m),gpuF,FH{tr(2,m)});            
            x=fftGPU(x,tr(1,m),gpuF,F{tr(1,m)});
            x=bsxfun(@times,x,et{2}{m});          
            x=ifftGPU(x,tr(1,m),gpuF,FH{tr(1,m)});
            if any(et{5}{m}(:)==1);x=bsxfun(@times,x,1-et{5}{m})+bsxfun(@times,flipping(x,et{4}{m}),et{5}{m});end%FLIPPING FOR LARGER THAN 90DEG ROTATIONS

        end 
    end
    if sumR;x=multDimSum(x,4:ndT);end 
end

x=resampling(x,N,2);
