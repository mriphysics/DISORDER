function sH=buildShearlet(N,J,gpu,K,DFtype,QMFtype,full,complex)

%BUILDSHEARLET computes a shearlet system. The code is based on [1] G 
%Kutyniok, et al, "ShearLab 3D: Faithful digital shearlet transforms based 
%on compactly supported shearlets," ACM Trans. Math. Softw., 
%42(1):5:1–5:42, 2016 
%   SH=BUILDSHEARLET(N,J,{GPU},{K},{DFTYPE},{QMFTYPE},{FULL},{COMPLEX})
%   * N are the spatial dimensions of the system
%   * J is the number of scales of the desired shearlet system. It has to 
%   be >= 1
%   * {GPU} determines whether to use gpu computations
%   * {K} is a 1xJ array, specifying the level of shearing occuring on each 
%   scale. Each entry of shearLevels has to be >= 0. A shear level of K 
%   means that the generating shearlet is sheared 2^K times in each 
%   direction for each cone (2D) / both types of shearing for all three
%   pyramids (3D). 2D example: If J=3 and K=[1 1 2], S will contain 
%   (2*(2*2^1+1))+(2*(2*2^1+1))+(2*(2*2^2+1))=38 shearlets. 3D example: If 
%   J=2 and K=[1 2], S will contain 3*((2*2^1)+1)^2+3*((2*2^2)+1)^2=318
%   shearlets. In both cases we omit the lowpass shearlet and translation). 
%   Note that it is recommended not to use the full shearlet system but to 
%   omit shearlets lying on the border of the second cone (2D) / second 
%   and third pyramid (3D) as they only slightly differ from shearlets on 
%   the border of the first cone (2) / first and second pyramid (3D). It
%   defaults to ceil((1:J)/2).
%   * {DFTYPE} is a 2D directional filter that serves as the basis of the 
%   directional 'component' of the shearlets. The default choice is 
%   'dmaxflat4'. For small sized inputs or very large systems, the default 
%   directional filter might be too large. In this case, it is recommended 
%   to use 'cd'. Options are 'haar': the "Haar" filters / 'vk': McClellan 
%   transformed of the filter from the VK book / 'ko': orthogonal filter 
%   in the Kovacevic's paper / 'kos': smooth 'ko' filter / 'lax': 17x17 
%   by Lu, Antoniou and Xu / 'sk': 9x9 by Shah and Kalker / 'cd': 7 and 
%   9 McClellan transformed by Cohen and Daubechies / 'dvmlp': regular 
%   linear phase biorthogonal filter with 3 dvm / 'sinc': ideal filter 
%   (no perfect recontruction!) / 'dmaxflat': diamond maxflat filters 
%   obtained from a three stage ladder
%   * {QMFTYPE} is a 1D quadrature mirror filter defining the 'wavelet'
%   component of the shearlets, the string should end up in an integer P
%   which relates to the support and vanishing moments of the wavelets.
%   Options are 'maxflat[0-9]': maximally flat symmetric 4+2P+1-tap low
%   pass filter (default  maxflat2) / 'Daubechies[0-9]': minimal phase 
%   filters that generate wavelets which have a minimal support for a 
%   given number of vanishing moments. The length is 2P+2 and the number of
%   vanishing moments is P+1 / 'Beylkin0' places roots for the frequency 
%   response function close to the Nyquist frequency on the real axis / 
%   'Coiflet[0-4]' are designed to give both the mother and father 
%   wavelets 2(P+1) vanishing moments / 'Symmlet[0-6]' are wavelets within 
%   a minimum size support for a given number of vanishing moments, but 
%   they are as symmetrical as possible, as opposed to the Daubechies 
%   filters which are highly asymmetrical. P+4 specifies the number of 
%   vanishing moments and the size of the support is 2(P+4) / 
%   'Vaidyanathan0' gives an exact reconstruction, but does not satisfy any 
%   moment condition. It is optimized for speech coding / 'Battle[0-2]' 
%   generates spline orthogonal wavelet basis. 2P+1 gives the degree of the 
%   spline. The number of vanishing moments is 2(P+1)
%   * {FULL} is a logical value that determines whether a full shearlet 
%   system is computed or if shearlets lying on the border of the second 
%   cone (2D) / second and third pyramid (3D) are omitted. The default and 
%   recommended value is 0.
%   * {COMPLEX} indicates whether complex shearlets are generated, so both
%   sides of the spectrum are not treated symmetrically. Defaults to 0
%   ** SH is a structure containing the specified shearlet system. It 
%   contains the following fields:
%       .S: A NxM array of M 2D shearlets in the frequency domain (starting
%   at frequency 0).
%       .K: The respective input argument is stored here.
%       .full: The respective input argument is stored here.
%       .idxs: A Mx3 (2D) / Mx4 (3D) array, specifying each shearlet in 
%   the system in the format [cone scale shearing] (2D) / [cone scale 
%   shearing1 shearing2] (3D). The vertical cone in the time domain is 
%   indexed by 1 while the horizontal cone is indexed by 2. Note that the 
%   values for scale and shearing are limited by specified number of scales 
%   and shear levels. The lowpass shearlet is indexed by [0 0 0].
%       .dfw: A N sized dual frame weights matrix containing the absolute 
%   and squared sum over all shearlets stored in sH.S. These weights are
%   needed to compute the dual frame during reconstruction.
%       .rms: A 1xM array containing the root mean squares (L2-norm divided 
%   by sqrt(prod(N)) of all shearlets stored in sH.S. These values can be 
%   used to normalize shearlet coefficients to make them comparable.

if nargin<3 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
if nargin<4 || isempty(K);K=ceil((1:J)/2);end
if nargin<5 || isempty(DFtype);DFtype='dmaxflat4';end
if nargin<6 || isempty(QMFtype);QMFtype='maxflat2';end
if nargin<7 || isempty(full);full=0;end
if nargin<8 || isempty(complex);complex=0;end%Experimental, not sure it works in 3D, the problem could be in getIndFromIdxs below

%INITIALIZE
ND=length(N);
cf=nchoosek(1:ND,2);
if ND==3;cf(3,:)=[3 2];end%It was the convention of previous code...
idxs=cell(1,ND-1);
for n=2:ND
    if n~=ND;idxs{n-1}=getShearletIdxs(n,K,1,complex);
    else idxs{n-1}=getShearletIdxs(n,K,full,complex);
    end
end

%WE BUILD THE DIRECTIONAL, QUADRATURE MIRROR (SCALING)
DF=buildDirectionalFilter(DFtype);
SF=buildQuadratureMirrorFilter(QMFtype);

%CONSTRUCT 2D SHEARLETS
for c=1:size(cf,1)
    F2D=prepareFilters(N(cf(c,:)),DF,SF);            
    [S{c},rms,dfw]=getShearlets(F2D,N(cf(c,:)),idxs{1});
    perm([cf(c,:) setdiff(1:3,cf(c,:)) 4])=[1:2 4 3];
    S{c}=permute(S{c},perm);
end

%COMBINE 2D SHEARLETS
if ND==3;[S{1},rms,dfw]=combineShearlets(S,idxs{2});end
sH=struct('S',S{1},'K',K,'full',full,'idxs',idxs{ND-1},'dfw',dfw,'rms',rms);

function F=prepareFilters(N,DF,SF)
    [DF,SF1,WF,SF2]=checkFilterFeasibility(N,K,DF,SF);
    if gpu;[DF,SF1,WF,SF2]=parUnaFun({DF,SF1,WF,SF2},@gpuArray);end
    for m=1:2        
        [F{m}.WD,F{m}.BP,F{m}.LP]=buildWedgeBandpassAndLowpassFilters(N,K,DF,SF1,WF,SF2,complex);  
        N=circshift(N,[0 -1]);
    end
end

function [S,rms,dfw]=getShearlets(F,N,idxs)
    NS=size(idxs,1);S=zeros([N NS],'like',F{1}.LP);
    for s=1:NS
        if idxs(s,1)==0;S(:,:,s)=F{1}.LP;
        elseif idxs(s,1)==1;S(:,:,s)=F{1}.WD{K(idxs(s,2))+1}(:,:,-idxs(s,4)+2^K(idxs(s,2))+1).*conj(F{1}.BP(:,:,idxs(s,2),idxs(s,3)));
        else S(:,:,s)=permute(F{2}.WD{K(idxs(s,2))+1}(:,:,idxs(s,4)+2^K(idxs(s,2))+1).*conj(F{2}.BP(:,:,idxs(s,2),idxs(s,3))),[2 1]);
        end
    end    
    [rms,dfw]=computeNormalization(S,N);
    rms=permute(rms,[1 2 4 3]);
end

function [S,rms,dfw]=combineShearlets(F,idxs)
    NS=size(idxs,1);S=zeros([N NS],'like',F{1});    
    for s=1:NS
        if idxs(s,1)==0;S=dynInd(S,s,4,bsxfun(@times,dynInd(F{1},size(F{1},4),4),dynInd(F{3},size(F{3},4),4)));
        elseif idxs(s,1)==1;S=dynInd(S,s,4,bsxfun(@times,dynInd(F{1},getIndFromIdxs(1,idxs(s,2),idxs(s,3),idxs(s,4)),4),dynInd(F{3},getIndFromIdxs(1,idxs(s,2),idxs(s,3),idxs(s,5)),4)));
        elseif idxs(s,1)==2;S=dynInd(S,s,4,bsxfun(@times,dynInd(F{2},getIndFromIdxs(1,idxs(s,2),idxs(s,3),idxs(s,4)),4),dynInd(F{3},getIndFromIdxs(2,idxs(s,2),idxs(s,3),idxs(s,5)),4)));
        else S=dynInd(S,s,4,bsxfun(@times,dynInd(F{1},getIndFromIdxs(2,idxs(s,2),idxs(s,3),idxs(s,4)),4),dynInd(F{2},getIndFromIdxs(2,idxs(s,2),idxs(s,3),idxs(s,5)),4)));
        end
    end
    [rms,dfw]=computeNormalization(S,N);
end

function [rms,dfw]=computeNormalization(S,N)
    NND=length(N);
    S=abs(S).^2;
    rms=sqrt(multDimSum(S,1:ND)/prod(N));
    dfw=sum(S,NND+1);
end

function ind=getIndFromIdxs(cone,scale,half,shearing)
    ind=(1+complex)*sum(2.^(K+1)+1)*(cone-1)+(2^K(scale))+half+shearing;
    if scale>1;ind=ind+(1+complex)*sum(2.^(K(1:(scale-1))+1)+1);end
end

end