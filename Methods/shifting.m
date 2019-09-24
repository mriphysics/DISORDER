function x=shifting(x,H,f,zp)

% SHIFTING applies a shift to an image in Fourier domain
%   X=SHIFTING(X,H,{F})
%   * X is the array to be shifted
%   * H is a cell array with the shifts to be applied along different 
%   dimensions
%   * {F} indicates if the data is in the Fourier domain already, defaults 
%   to 0
%   * {ZP} serves to zero pad the output to emulate non-circular shifts, 
%   only implemented for circular shifts at the moment
%   ** X is the shifted image
%

if nargin<3 || isempty(f);f=0;end
if nargin<4 || isempty(zp);zp=0;end

gpu=isa(x,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

NH=length(H);
NX=size(x);NX(end+1:NH)=1;

rea=isreal(x);
ND=length(NX);
for m=1:NH
    HH=H{m};
    NXX=ones(1,ND);NXX(m)=NX(m);
    if ~isempty(HH)
        if numel(HH)==1 && all(mod(HH(:),1))==0
            if NXX(m)~=1
                rGrid=generateGrid(NXX(m),0,NXX(m),0);
                indGrid=circshift(rGrid{1},HH);            
                x=dynInd(x,indGrid,m);
                if zp
                    if indGrid(1)<NXX(m)/2;x=dynInd(x,find(indGrid==1):NXX(m),m,0);else x=dynInd(x,1:find(indGrid==NXX(m)),m,0);end
                end
            end
        else            
            kGrid=generateGrid(NXX,gpu,2*pi,ceil((NXX+1)/2));
            kGrid=-1i*ifftshift(kGrid{m},m);
            HH=exp(bsxfun(@times,HH,kGrid));

            if ~f;x=fftGPU(x,m,gpuF);end%Otherwise we asume it is already in Fourier domain
            x=bsxfun(@times,x,HH);
            if ~f;x=ifftGPU(x,m,gpuF);end
            if rea;x=real(x);end        
        end
    end
end
