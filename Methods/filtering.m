function x=filtering(x,H,c,f,quick,F,FH)

% FILTERING filters an array in Fourier domain
%   X=FILTERING(X,H,{C},{F})
%   * X is the array to be filtered
%   * H is the filter to be applied
%   * {C} indicates if the filter is to be applied in the Cosine domain 
%   (for Neumann instead of periodic boundary conditions), defaults to 0, 
%   i.e., the filter is applied in the Fourier domain
%   * {F} indicates if the data is in the Fourier domain already, defaults 
%   to 0
%   * {QUICK} serves to launch quick resampling when the data is real, with
%   1 the filter hast to be made yet compatible with real data, with 2 it
%   is already compatible
%   * {F} is a cell containing the unitary DFT matrix along different 
%   dimensions
%   * {FH} is a cell containing the unitary IDFT matrix along different 
%   dimensions
%   ** X is the filtered image
%

if nargin<3 || isempty(c);c=0;end
if nargin<4 || isempty(f);f=0;end
if nargin<5 || isempty(quick);quick=0;end

gpu=isa(x,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

NH=size(H);
nDimsH=ndims(H);
comp=~isreal(x) | ~isreal(H) | ~quick;%NOTE THE CODE MAY NOT WORK FOR COMP=1

if length(c)==1;c=c*ones(1,nDimsH);end

if gpuF~=2 || comp;quick=0;end%For CPU data matrix multiplication is not very efficient probably

if nargin<7
    if quick;[F,FH]=buildStandardDFTM(NH,0,gpu,~comp);else F=cell(1,nDimsH);FH=cell(1,nDimsH);end
end

for m=1:nDimsH
    if NH(m)~=1 && ~f
        if ~c(m);x=fftGPU(x,m,gpuF,F{m});
        else x=fctGPU(x,m,gpuF);
        end
    end
end

if quick==1 && ~comp%This makes the filter compatible with real data
    for m=1:nDimsH
        if ~c(m)            
            Nre=ceil((NH(m)+1)/2);
            ind=repmat(1:Nre,[2 1]);
            ind(2)=[];
            if mod(NH(m),2)==0;ind(end)=[];end
            H=dynInd(H,ind(:),m);
        end
    end
end
x=bsxfun(@times,x,H);

for m=1:nDimsH
    if NH(m)~=1 && ~f
        if ~c(m);x=ifftGPU(x,m,gpuF,FH{m});
        else x=ifctGPU(x,m,gpuF);
        end
    end
end