function x=resampling(x,Nres,fo,mirror,tr,quick)

% RESAMPLING resamples a given array using the FFT
%   X=RESAMPLING(X,NRES,{FO},{MIRROR})
%   * X is the array to be resampled
%   * NRES is the new grid size
%   * {FO} determines whether the input/output is in k-space (1), shifted
%   k-space (2) or space (0), defaults to 0
%   * {MIRROR} determines whether to mirror the image along a given
%   dimension, defaults to 0 for all dimensions 
%   * {QUICK} serves to launch quick resampling when the data is real
%   * {TR} trims the singleton dimensions on Nres 
%   ** X is the resampled array
%

if nargin<3 || isempty(fo);fo=0;end
if nargin<5 || isempty(tr);tr=1;end
if nargin<6 || isempty(quick);quick=0;end

BlSz=1e6;

gpu=isa(x,'gpuArray');if gpu;gpuF=2;else gpuF=0;end
comp=~isreal(x);

N=size(x);nDimsIn=length(N);
if tr==0;indTrim=[];else indTrim=find(Nres~=1,1,'last');end
if ~isempty(indTrim) && indTrim~=1;Nres(indTrim+1:end)=[];end
nDimsOu=length(Nres);

if nargin<4 || isempty(mirror);mirror=single(zeros(1,nDimsOu));end

if tr==1;assert(nDimsOu<=nDimsIn,'Resampling dimensionality (%d) is larger than image dimensionality (%d)',nDimsOu,nDimsIn);elseif nDimsIn<nDimsOu;N(end+1:nDimsOu)=1;end
Nor=N(1:nDimsOu);mirror=mirror(1:nDimsOu);mirrorin=mirror;mirrorin(:)=0;

NorM=Nor+(mirror==1).*Nor;NresM=Nres+(mirror==1).*Nres;
Nmin=min(NorM,NresM);Nmax=max(NorM,NresM);

zeroF=ceil((Nmax+1)/2);
orig=zeroF-ceil((Nmin-1)/2);
fina=zeroF+floor((Nmin-1)/2);
orig(mirror==2)=1;
fina(mirror==2)=Nmin(mirror==2);

if gpuF~=2 || comp;quick=0;end%For CPU data matrix multiplication is not very efficient probably, for not complex data, not completely clear how to operate yet
for m=1:nDimsOu
    if Nor(m)~=Nres(m)
        NNres=[Nres(1:m-1) NresM(m) N(m+1:end)];          
        mirrorin(m)=mirror(m);
        if mirror(m)~=2;x=mirroring(x,mirrorin==1,1);end%MIRRORING FUNCTION NOT INCLUDED
        if ~fo
            if mirror(m)==2
                x=fctGPU(x,m,gpuF)/sqrt(NorM(m));%FCT FUNCTION NOT INCLUDED
            elseif ~quick
                x=fftGPU(x,m,gpuF)/NorM(m);
            else               
                F=build1DFTM(NorM(m),0,gpu,~comp)/NorM(m);
                x=fftGPU(x,m,gpuF,F,~comp);
            end
        end
        if ~quick || fo~=0 || mirror(m)==2 || comp
            if fo~=2 && mirror(m)~=2
                if numel(x)<1e8
                    x=fftshift(x,m);
                else
                    perm=1:nDimsIn;perm([1 m])=[m 1];
                    x=permute(x,perm);
                    NX=size(x);
                    x=x(:,:);
                    N2=size(x,2);
                    for s=1:BlSz:N2;vS=s:min(s+BlSz-1,N2);x(:,vS)=fftshift(x(:,vS),1);end                  
                    x=reshape(x,NX);
                    x=permute(x,perm);
                end 
            end
            xRes=zeros(NNres,'like',x);
            if Nor(m)<Nres(m);x=dynInd(xRes,orig(m):fina(m),m,x);else x=dynInd(x,orig(m):fina(m),m);end
            if fo~=2 && mirror(m)~=2;x=ifftshift(x,m);end           
            if ~fo
                if mirror(m)==2;x=ifctGPU(x,m,gpuF)*sqrt(NresM(m));else x=ifftGPU(x,m,gpuF)*NresM(m);end                 
            end
        else
            [~,FH]=build1DFTM(NresM(m),0,gpu,~comp);
            if NresM(m)>NorM(m)
                FH=FH(:,1:NorM(m))*NresM(m);
                if mod(NorM(m),2)==0;FH(:,end)=FH(:,end)/2;end
            else
                x=dynInd(x,1:NresM(m),m)*NresM(m);
            end
            x=ifftGPU(x,m,gpuF,FH,~comp);
        end
        if mirror(m)~=2;x=mirroring(x,mirrorin==1,0);end
        mirrorin(m)=0;
    end
end
if ~comp;x=real(x);end
