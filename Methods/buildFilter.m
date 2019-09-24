function x=buildFilter(N,typ,sp,gpu,gibbsRing,c)

% BUILDFILTER builds a filter in the frequency domain
%   X=BUILDFILTER(N,TYP,{SP},{GPU},{GIBBSRING},{C})
%   * N are the dimensions of the filter
%   * TYP is the filter type
%   * {SP} is the physical spacing
%   * {GPU} determines whether to use gpu computations 
%   * {GIBBSRING} is a filter parameter to prevent Gibbs ringing, defaults 
%   to 0
%   * {C} indicates if the filter is to be applied in the Cosine domain 
%   (for Neumann instead of periodic boundary conditions), defaults to 0, 
%   i.e., the filter is applied in the Fourier domain
%   ** X is the filter profile
%

if nargin<3 || isempty(sp);sp=1;end
if nargin<4 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
if nargin<5 || isempty(gibbsRing);gibbsRing=0;end
if nargin<6 || isempty(c);c=0;end

if gpu;gpuF=2;else gpuF=0;end

ND=length(N);
if length(c)==1;c=c*ones(1,ND);end
if length(sp)==1;sp=sp*ones(1,ND);end
x=0;
if gpu;x=gpuArray(x);end
%%NOT SURE THAT THE SPACING IN KGRID IS EQUIVALENT TO THE SPACING IN
%%RGRID...
kGrid=generateGrid(N,gpu,2*1i*pi./sp,ceil((N+1)/2));
rGrid=generateGrid(N,gpu);for n=1:ND;rGrid{n}(:)=0;end

if strcmp(typ,'1stFinite')%First order continuous finite differences   
    for n=1:ND;x=bsxfun(@plus,x,kGrid{n});end   
elseif strcmp(typ,'2ndFinite')%Second order continuous finite differences
    for n=1:ND;x=bsxfun(@plus,x,kGrid{n}.^2);end
elseif strcmp(typ,'1stFiniteDiscreteForward')%First order forward finite differences
    assert(all(~c),'1st Order Finite Discrete Difference filter not defined for DCT\n');
    for n=1:ND
        if N(n)>1;rGrid{n}(1)=1;rGrid{n}(end)=-1;rGrid{n}=-rGrid{n}/(sp(n));end
    end
elseif strcmp(typ,'1stFiniteDiscreteBackward')%First order forward finite differences
    assert(all(~c),'1st Order Finite Discrete Difference filter not defined for DCT\n');
    for n=1:ND
        if N(n)>1;rGrid{n}(1)=-1;rGrid{n}(2)=1;rGrid{n}=-rGrid{n}/(sp(n));end
    end
elseif strcmp(typ,'2ndFiniteDiscrete')%This corresponds, for instance, to 
    %the cosine profile in DC Ghiglia and LA Romero, "Robust 
    %two-dimensional weighted and unweighted phase unwrapping that uses
    %fast transforms and iterative methods", J Opt Soc Am A, 11(1):107-117, 
    %Jan 1994
    for n=1:ND
        if N(n)>2;rGrid{n}(1)=-2;rGrid{n}(2)=1;rGrid{n}(end)=1;rGrid{n}=rGrid{n}/(sp(n)^2);end
    end
elseif strcmp(typ,'FractionalFiniteDiscrete')
    for n=1:ND
        if N(n)>1;x=bsxfun(@plus,x,(2*abs(sin(imag(kGrid{n})/2))).^gibbsRing);end        
    end
elseif strcmp(typ,'FractionalFiniteDiscreteIso') || strcmp(typ,'FractionalFiniteDiscreteIsoNorm') || strcmp(typ,'FractionalIso')
    for n=1:ND
        if N(n)>1;x=bsxfun(@plus,x,abs(kGrid{n}).^2);end
    end
    x=sqrt(x);
    if strcmp(typ,'FractionalFiniteDiscreteIso') || strcmp(typ,'FractionalFiniteDiscreteIsoNorm')
        x(x>pi)=pi;
        if strcmp(typ,'FractionalFiniteDiscreteIso');mult=2;else mult=1;end
        x=(mult*abs(sin(x/2)));
    end
    x=x.^gibbsRing;
elseif strcmp(typ,'tukeyIso')
    kk=single(zeros(prod(N),ND));
    if gpu>0;kk=gpuArray(kk);end   
    for m=1:ND
        rep=N;rep(m)=1;
        kkAux=repmat(kGrid{m},rep);
        kk(:,m)=kkAux(:);
    end;kkAux=[];   
    kk=reshape(kk,[N ND]);  
    dimkk=ndims(kk);   
    %kkrad=sqrt(sum(abs(kk.*conj(kk)),dimkk));%Radial k-coordinates
    %To save memory we perform previous instruction as separate steps
    kk=abs(kk); 
    kk=kk.^2;  
    kk=sum(kk,dimkk);   
    kk=sqrt(kk);
    N(end+1:2)=1;x=single(ones(N));
    if gpu>0;x=gpuArray(x);end
    alpha=1-gibbsRing;
    if gibbsRing~=0
        fkk=0.5*(1+cos(pi*((kk-pi*alpha)/((1-alpha)*pi))));
        x(kk>=pi*alpha)=fkk(kk>=pi*alpha);
    end
    x(kk>=pi)=0;kk=[];
elseif strcmp(typ,'tukey')
    N(end+1:2)=1;x=single(ones(N));
    if gpu>0;x=gpuArray(x);end
    gibbsRing(end+1:ND)=gibbsRing(end);
    for m=1:ND
        Naux=ones(1,ND);Naux(m)=N(m);spaux=ones(1,ND);spaux(m)=sp(m);        
        x=bsxfun(@times,x,buildFilter(Naux,'tukeyIso',spaux,gpu,gibbsRing(m)));
    end
    x=fftshift(x);
elseif strcmp(typ,'CubicBSpline')
    x=1;
    for n=1:ND
        x=bsxfun(@times,x,6*bsxfun(@rdivide,sinc(imag(kGrid{n}/(2*pi))).^4,4+2*cos(imag(kGrid{n}))));
    end
else
    error('Unknown filter %s',typ);
end


if strcmp(typ,'1stFiniteDiscreteForward') || strcmp(typ,'1stFiniteDiscreteBackward') || strcmp(typ,'2ndFiniteDiscrete')
    for n=1:ND
        if numel(rGrid{n})>1
            rGrid{n}=fftGPU(rGrid{n},n,gpuF);
            if c(n);rGrid{n}=real(rGrid{n});end
            rGrid{n}=fftshift(rGrid{n},n);
            x=bsxfun(@plus,x,rGrid{n});
        end
    end
end
x=ifftshift(x);
if any(c)%We assume the filter is symmetric
    N=size(x);    
    assert(isreal(x) || any(~c),'The filter is not real, so it cannot be applied in the cosine domain');
    assert(all(mod(N(c==1),2)==0 | N(c==1)==1),'The filter does not have an even size in all dimensions as its size is%s',sprintf('%d',N)); 
    v=cell(1,ND);
    for m=1:ND
        if c(m);v{m}=1:ceil(N(m)/2);else v{m}=1:N(m);end
    end
    x=dynInd(x,v,1:ND);
end

if strcmp(typ,'FractionalFiniteDiscreteIso') || strcmp(typ,'FractionalFiniteDiscreteIsoNorm') || strcmp(typ,'FractionalIso') && gibbsRing<0;x(1)=max(x(2:end));end%To prevent numerical instabilities
