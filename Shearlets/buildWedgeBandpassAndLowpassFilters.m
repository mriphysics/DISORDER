function [WD,BP,LP]=buildWedgeBandpassAndLowpassFilters(N,K,DF,SF1,WF,SF2,complex)

%BUILDWEDGEBANDPASSANDLOWPASSFILTERS builds the filters that comprise the
%shearlet system
%   [WD,BP,LP]=BUILDWEDGEBANDPASSANDLOWPASSFILTERS(N,K,DF,SF1,WF,SF2,{COMPLEX}) 
%   builds the filters that comprise the shearlet system. Page and equation 
%   numbers refer to [1] G Kutyniok, et al, "ShearLab 3D: Faithful Digital 
%   Shearlet Transforms Based on Compactly Supported Shearlets," ACM TMS 
%   42,5, 2016
%   * N are the dimensions of the shearlet
%   * K specifies the levels of shearing for each scale
%   * DF is the directional filter
%   * SF1 is the first scaling filter
%   * WF is the wavelet filter
%   * SF2 is the second scaling filter
%   * {COMPLEX} indicates whether complex shearlets are generated, so both
%   sides of the spectrum are not treated symmetrically. Defaults to 0
%   ** WD is the wedge filter
%   ** BP is the bandpass filter
%   ** LP is the lowpass filter
%

if nargin<7;complex=0;end

gpu=isa(DF,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

%INITIALIZATION
[DF,SF1,WF,SF2]=parUnaFun({DF,SF1,WF,SF2},@single);
N(end+1:2)=1;
J=length(K);%Number of scales
BP=zeros([N J],'like',DF);%These filters partition the frequency plane into different scales
Kmax=max(K)+1;
WD=cell(1,Kmax);%These filters partition the frequency plane into different directions
DF=DF/sum(abs(DF(:)));%This filter corresponds to the time domain representation of the trigonometric polynomial P (see equation (13) on page 12)

%1D HIGH AND LOWPASS FILTERS AT DIFFERENT SCALES
FHP=cell(1,J);%High pass filter, FHP{J}=g_1 and FHP{1}=g_J (compare page 11) 
FLP1=cell(1,J);%Low pass filter 1, FLP1{J}=h_1 and FLP1{1}=h_J (compare page 11)
FLP2=cell(1,Kmax);%Low pass filter 2, typically we have FLP2{max(K)+1}=FLP1{J}, i.e. FLP2{J}=h_1 (compare page 11)

%INITIALIZE WAVELET HIGHPASS AND LOWPASS FILTERS
FHP{J}=WF;%This filter is typically chosen to form a quadrature mirror filter pair with SF1 and corresponds to g_1 on page 11
FLP1{J}=SF1;%This filter corresponds to h_1 on page 11
FLP2{Kmax}=SF2;%This filter is typically chosen to be equal to SF1 and provides the y-direction for the tensor product constructing the 2D wavelet filter w_j on page 14

%COMPUTE WAVELET HIGH- AND LOWPASS FILTERS ASSOCIATED WITH A 1D DIGITAL WAVELET TRANSFORM ON J SCALES, E.G. WE COMPUTE h_1 TO h_J AND g_1 TO g_J (COMPARE PAG 11)
for j=J-1:-1:1
    FLP1{j}=conv(FLP1{J},upsampling(FLP1{j+1},2));
    FHP{j}=conv(FLP1{J},upsampling(FHP{j+1},2));
end
for k=Kmax-1:-1:1
    FLP2{k}=conv(FLP2{Kmax},upsampling(FLP2{k+1},2));
end

%CONSTRUCT BANDPASS FILTERS FOR SCALES 1 TO J
for j=1:J;BP(:,:,j)=resampling(FHP{j},N,2);end
BP=fft2D(BP);

if complex
    %HALF FILTER
    HF=zeros(N,'like',BP);
    HF(:,1:ceil(N(2)/2))=1;
    HF(:,:,1,2)=1-HF;
    HF=ifftshift(HF,2);
    BP=bsxfun(@times,BP,HF);
end

%CONSTRUCT WEDGE FILTERS FOR ACHIEVING DIRECTIONAL SELECTIVITY. AS THE
%ENTRIES IN K DESCRIBE THE NUMBER OF DIFFERENTLY SHARED ATOMS AT A CERTAIN
%SCALE, A DIFFERENT SET OF WEDGE FILTERS HAS TO BE CONSTRUCTED FOR EACH
%VALUE IN K
for k=unique(K)
    %PREALLOCATE A TOTAL OF FLOOR(2^(k+1)+1) WEDGE FILTERS, THE NUMBER OF
    %DIFFERENT DIRECTIONS OF SHEARLET ATOMS ASSOCIATED WITH THE HORIZONTAL
    %/ VERTICAL FREQUENCY CONES    
    WD{k+1}=zeros([N floor(2^(k+1)+1)],'like',DF);%Plus one for one unsheared shearlet
    
    %UPSAMPLE DIRECTIONAL FILTER IN THE VERTICAL DIRECTION. BY UPSAMPLING
    %THE DIRECTIONAL FILTER IN THE TIME DOMAIN, WE CONSTRUCT REPEATING
    %WEDGES IN THE FREQUENCY DOMAIN (COMPARE abs(fftshift(fft2(ifftshift(DFUPS))))
    %AND abs(fftshift(fft2(ifftshift(DF)))))
    DFups=upsampling(DF,1,2^(k+1));
    
    %REMOVE HIGH FREQUENCIES ALONG THE VERTICAL DIRECTION IN THE FREQUENCY
    %DOMAIN BY CONVOLVING THE UPSAMPLED DIRECTIONAL FILTER WITH A LOWPASS
    %FILTER. WE REMOVE ALL BUT THE CENTRAL WEDGE IN THE FREQUENCY DOMAIN
    %WDAUX CORRESPONDS TO conv(p_j,h_(J-j*alpha_j/2)') IN THE LANGUAGE OF 
    %THE PAPER. TO SEE THIS CONSIDER THE DEFINITION OF p_j ON PAGE 14, THE 
    %DEFINITION OF w_j ON THE SAME PAGE AND THE DEFINITION OF THE DIGITAL 
    %SHEARLET FILTER ON PAG 15. FURTHERMORE, THE g_j PART OF THE 2D 
    %WAVELET FILTER w_j IS INVARIANT TO SHEARINGS, HENCE IT SUFFICES TO 
    %APPLY THE DIGITAL SHEAR OPERATOR TO WDAUX (COMPARE EQUATION (22)). 
    %WE UPSAMPLE WEDGE FILTER IN X-DIRECTION. THIS OPERATION CORRESPONDS TO
    %THE UPSAMPLING EQ (21) ON PAG 15
    WDups=upsampling(resampling(conv2(DFups,FLP2{Kmax-k}'),N,2),2,2^k);  
    
    %CONVOLVE WEDGE FILTER WITH LOWPASS FILTER, AGAIN FOLLOWING EQ (21) ON 
    %PAG 14.
    FLPaux=resampling(FLP2{Kmax-max(k-1,0)},size(WDups),2);    
    if k>=1;WDups=ifft2D(fft2D(WDups).*fft2D(FLPaux));end
    FLPauxfl=flip(FLPaux,2);    
  
    %TRAVERSE ALL DIRECTIONS OF THE UPPER PART OF THE LEFT HORIZONTAL
    %FREQUENCY CONE
    for l=-2^k:2^k
        %RESAMPLE WDUPS AS GIVEN IN EQ (22) ON PAG 15
        WDupssh=shearing(WDups,l,[2 1]);      
        %CONVOLVE AGAIN WITH FLIPPED LPF, AS REQUIRED BY EQ (22) ON PAG 15
        if k>=1;WDupssh=ifft2D(fft2D(WDupssh).*fft2D(FLPauxfl));end
        %OBTAIN DOWNSAMPLED AND RENORMALIZED AND SHEARED WEDGE FILTER IN 
        %THE FREQUENCY DOMAIN, ACCORDING TO EQ (22) ON PAG 15
        WDupssh=(2^k)*WDupssh(:,1:2^k:(2^k)*N(2));
        WDupssh=fft2D(WDupssh);
        WD{k+1}(:,:,2^k+1-l)=WDupssh;
    end
end

%COMPUTE LOW PASS FILTER OF SHEARLET SYSTEM
LP=fft2D(resampling(FLP1{1}'*FLP1{1},N,2));    

function x=fft2D(x)
    for m=1:2;x=ifftshift(x,m);x=fftGPU(x,m,gpuF);end
end

function x=ifft2D(x)
    for m=1:2;x=ifftGPU(x,m,gpuF);x=fftshift(x,m);end
end

end
