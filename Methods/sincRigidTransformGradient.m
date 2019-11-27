function [G,GB,GC]=sincRigidTransformGradient(xB,et,etg,F,FH)

%SINCRIGIDTRANSFORMGRADIENT obtains the gradient of the sinc rigid 
%transform of the volumes
%   G=SINCRIGIDTRANSFORMGRADIENT(XB,ET,ETG,{F},{FH})
%   * XB are the volumes before the first, second and third rotations, 
%   before the translation and original (for the padding dimensions)
%   * ET are the transform factors
%   * ETG are the transform gradient factors
%   * {F} contains discrete Fourier transform matrices
%   * {FH} contains inverse discrete Fourier transform matrices
%   ** G is the gradient of the transformed image
%   ** GB is the gradient of the first, second and third rotations before 
%   applying the translation
%   ** GC is the gradient of the first, first and second rotations before
%   applying the second, third and third rotations respectively
%

if nargin<4 || isempty(F);F={[],[],[]};end
if nargin<5 || isempty(FH);FH={[],[],[]};end

G=cell(1,6);
if nargout>=2;GB=cell(1,3);end
if nargout>=3;GC=cell(1,3);end

N=size(xB{5});N(end+1:3)=1;N=N(1:3);

%TRANSLATION PARAMETERS
for m=1:3  
    x{1}=bsxfun(@times,xB{4},etg{1}{m});
    if iscell(et{1})  
        for n=1:3;x{1}=bsxfun(@times,x{1},et{1}{n});end
    end
    for n=1:3;x{1}=ifftGPU(x{1},n,FH{n});end    
    G{m}=x{1};
end

%FIRST ROTATION
x{1}=bsxfun(@times,xB{1},et{2}{1});
x{2}=bsxfun(@times,xB{1},etg{2}{1});
for m=1:2
    x{m}=ifftGPU(x{m},1,FH{1});
    x{m}=fftGPU(x{m},2,F{2});
end
x{2}=bsxfun(@times,x{1},etg{3}{1})+bsxfun(@times,x{2},et{3}{1});
x{1}=bsxfun(@times,x{1},et{3}{1});    
for m=1:2
    x{m}=ifftGPU(x{m},2,FH{2});       
    x{m}=fftGPU(x{m},1,F{1});
end
x{1}=bsxfun(@times,x{1},etg{2}{1})+bsxfun(@times,x{2},et{2}{1});
x{1}=ifftGPU(x{1},1,FH{1});

if any(et{5}{2}(:)==1);x{1}=bsxfun(@times,x{1},1-et{5}{2})+bsxfun(@times,flipping(x{1},et{4}{2}),et{5}{2});end
x{1}=fftGPU(x{1},3,F{3});
if nargout>=3;GC{1}=x{1};end
x{1}=bsxfun(@times,x{1},et{2}{2});
x{1}=ifftGPU(x{1},3,FH{3});
x{1}=fftGPU(x{1},1,F{1});
x{1}=bsxfun(@times,x{1},et{3}{2});
x{1}=ifftGPU(x{1},1,FH{1});
x{1}=fftGPU(x{1},3,F{3});
x{1}=bsxfun(@times,x{1},et{2}{2});
x{1}=ifftGPU(x{1},3,FH{3});

if any(et{5}{3}(:)==1);x{1}=bsxfun(@times,x{1},1-et{5}{3})+bsxfun(@times,flipping(x{1},et{4}{3}),et{5}{3});end
x{1}=fftGPU(x{1},2,F{2});
if nargout>=3;GC{2}=x{1};end
x{1}=bsxfun(@times,x{1},et{2}{3});
x{1}=ifftGPU(x{1},2,FH{2});
x{1}=fftGPU(x{1},3,F{3});
x{1}=bsxfun(@times,x{1},et{3}{3});
x{1}=ifftGPU(x{1},3,FH{3});
x{1}=fftGPU(x{1},2,F{2});
x{1}=bsxfun(@times,x{1},et{2}{3});

for m=1:2:3;x{1}=fftGPU(x{1},m,F{m});end
if nargout>=2;GB{1}=x{1};end
if ~iscell(et{1})
    x{1}=bsxfun(@times,x{1},et{1});
else
    for m=1:3;x{1}=bsxfun(@times,x{1},et{1}{m});end
end
for m=3:-1:1;x{1}=ifftGPU(x{1},m,FH{m});end
G{4}=x{1};

%SECOND ROTATION
x{1}=bsxfun(@times,xB{2},et{2}{2});
x{2}=bsxfun(@times,xB{2},etg{2}{2});
for m=1:2
    x{m}=ifftGPU(x{m},3,FH{3});
    x{m}=fftGPU(x{m},1,F{1});
end
x{2}=bsxfun(@times,x{1},etg{3}{2})+bsxfun(@times,x{2},et{3}{2});
x{1}=bsxfun(@times,x{1},et{3}{2});
for m=1:2
    x{m}=ifftGPU(x{m},1,FH{1});
    x{m}=fftGPU(x{m},3,F{3});
end
x{1}=bsxfun(@times,x{1},etg{2}{2})+bsxfun(@times,x{2},et{2}{2});
x{1}=ifftGPU(x{1},3,FH{3});

if any(et{5}{3}(:)==1);x{1}=bsxfun(@times,x{1},1-et{5}{3})+bsxfun(@times,flipping(x{1},et{4}{3}),et{5}{3});end
x{1}=fftGPU(x{1},2,F{2});
if nargout>=3;GC{3}=x{1};end
x{1}=bsxfun(@times,x{1},et{2}{3});
x{1}=ifftGPU(x{1},2,FH{2});
x{1}=fftGPU(x{1},3,F{3});
x{1}=bsxfun(@times,x{1},et{3}{3});
x{1}=ifftGPU(x{1},3,FH{3});
x{1}=fftGPU(x{1},2,F{2});
x{1}=bsxfun(@times,x{1},et{2}{3});

for m=1:2:3;x{1}=fftGPU(x{1},m,F{m});end
if nargout>=2;GB{2}=x{1};end
if ~iscell(et{1})
    x{1}=bsxfun(@times,x{1},et{1});
else
    for m=1:3;x{1}=bsxfun(@times,x{1},et{1}{m});end
end
for m=3:-1:1;x{1}=ifftGPU(x{1},m,FH{m});end
G{5}=x{1};

%THIRD ROTATION
x{1}=bsxfun(@times,xB{3},et{2}{3});
x{2}=bsxfun(@times,xB{3},etg{2}{3});
for m=1:2
    x{m}=ifftGPU(x{m},2,FH{2});
    x{m}=fftGPU(x{m},3,F{3});
end
x{2}=bsxfun(@times,x{1},etg{3}{3})+bsxfun(@times,x{2},et{3}{3});
x{1}=bsxfun(@times,x{1},et{3}{3});
for m=1:2
    x{m}=ifftGPU(x{m},3,FH{3});
    x{m}=fftGPU(x{m},2,F{2});
end
x{1}=bsxfun(@times,x{1},etg{2}{3})+bsxfun(@times,x{2},et{2}{3});

for m=1:2:3;x{1}=fftGPU(x{1},m,F{m});end
if nargout>=2;GB{3}=x{1};end
if ~iscell(et{1})
    x{1}=bsxfun(@times,x{1},et{1});
else
    for m=1:3;x{1}=bsxfun(@times,x{1},et{1}{m});end
end
for m=3:-1:1;x{1}=ifftGPU(x{1},m,FH{m});end
G{6}=x{1};

for n=1:length(G);G{m}=resampling(G{m},N,2);end
    
