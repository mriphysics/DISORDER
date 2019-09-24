function xou=regularize(x,R,pr)

%REGULARIZE   Applies a given regularization operator to the data
%   XOU=REGULARIZE(X,R,PR)
%   * X is the data before regularization
%   * R is the structure for regularization
%   * PR indicates whether to compute a prior-type regularizer
%   ** XOU is the data after regularization
%

if nargin<3;pr=[1 2];end
gpu=isa(x,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end

xou=x;xou(:)=0;
if ismember(1,pr)
    y=cell(1,5);
    if isfield(R,'Fd')
        y{1}=R.Fd.la*finiteDifferenceRegularizer(x,R.Fd.n,R.Fd.or,R.Fd.per,R.Fd.sp,R.Fd.sep);%All together although it does not give exactly the same
        %results as below
        %y{1}=finiteDifference(x,R.Fd.n,R.Fd.or,0,R.Fd.per,R.Fd.sp);%Forward
        %y{1}=R.Fd.la*finiteDifference(y{1},R.Fd.n,R.Fd.or,1,R.Fd.per,R.Fd.sp);%Backward    
        %y{1}=R.Fd.la*finiteDifferenceOld(x,R.Fd.n,R.Fd.or,R.Fd.per,R.Fd.sp,R.Fd.sep,gpu);%Old
    end%Finite differences
    if isfield(R,'Ti');y{2}=bsxfun(@times,R.Ti.la,x);end%Tikhonov
    if isfield(R,'Fo');y{3}=bsxfun(@times,R.Fo.la,filtering(x,conj(R.Fo.Fi).*R.Fo.Fi,R.Fo.mi));end
    if isfield(R,'Wa');y{4}=(R.Wa.la*R.Wa.p/2)*waveletTransform(bsxfun(@times,R.Wa.We,waveletTransform(x,R.Wa.wA)),R.Wa.wA,0);end%Wavelet
    if isfield(R,'Sh')
        xF=x;
        ND=numDims(x);
        for m=1:ND;xF=fftGPU(xF,m,gpuF);end
        NW=size(R.Sh.We,4);
        y{5}=xF;y{5}(:)=0;
        for w=1:NW
            sh=R.Sh.sH.S(:,:,:,w);
            we=R.Sh.We(:,:,:,w);
            if gpu;sh=gpuArray(sh);we=gpuArray(we);end
            xFF=bsxfun(@times,xF,conj(sh));
            for m=1:ND;xFF=ifftGPU(xFF,m,gpuF);end
            xFF=bsxfun(@times,xFF,we);
            for m=1:ND;xFF=fftGPU(xFF,m,gpuF);end
            xFF=bsxfun(@times,xFF,sh);
            y{5}=y{5}+xFF;
        end
        for m=1:ND;y{5}=ifftGPU(y{5},m,gpuF);end
        y{5}=(R.Sh.la*R.Sh.p/2)*y{5};        
        %y{5}=(R.Sh.la*R.Sh.p/2)*shearletTransform(bsxfun(@times,R.Sh.We,shearletTransform(x,R.Sh.sH)),R.Sh.sH,-1);
    end%Shearlet   
    for n=1:length(y)
        if ~isempty(y{n});xou=xou+y{n};end
    end
end
if ismember(2,pr)    
    y=cell(1,1);
    if isfield(R,'Tp');y{1}=bsxfun(@times,R.Tp.la,x);end%Tikhonov prior
    for n=1:length(y)
        if ~isempty(y{n});xou=xou+y{1};end
    end
end
if ismember(3,pr)
    y=cell(1,1);
    if isfield(R,'Tp');y{1}=bsxfun(@times,R.Tp.la,R.Tp.x0);end%Tikhonov prior
    for n=1:length(y)
        if ~isempty(y{n});xou=xou+y{1};end
    end
end
if ismember(4,pr)
    y=cell(1,1);
    if isfield(R,'Tp') && isfield(R.Tp,'x0')
        y{1}=(x-R.Tp.x0);
        y{1}=bsxfun(@times,R.Tp.la,abs(y{1}).^2);
    end     
    for n=1:length(y)
        if ~isempty(y{n});xou=xou+y{1};end
    end
end