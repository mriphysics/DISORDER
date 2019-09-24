function xou=decode(x,EH,E)

%DECODE   Applies a given decoding to the data
%   XOU=DECODE(X,EH,{E})
%   * X is the data before decoding
%   * EH is the decoding structure
%   * {E} is the encoding structure (memory demanding or replicated 
%   information may be stored in this structure)
%   ** XOU is the data after decoding
%

if isempty(EH);xou=x;return;end
if nargin<3;E=[];end
inBlSz=1e3;

if ~isfield(E,'bS')%Block sizes by default
    E.bS=[size(x,5) size(x,4)];E.dS=[size(x,5) size(x,4)];E.oS=[1 1];
    if isfield(E,'Sf');E.bS(2)=size(E.Sf,4);E.dS(2)=size(E.Sf,4);end
    if isfield(E,'Tr');E.bS(1)=E.NSe;E.dS(1)=E.NSe;end
end
if isfield(E,'Zf') && isempty(E.Zf) && isfield(E,'Sf');E.oS(2)=1;E.bS(2)=size(E.Sf,4);E.dS(2)=size(E.Sf,4);end

if isfield(EH,'Mb');x=bsxfun(@times,x,EH.Mb);end%Mask / Weighting
xou=x;
for a=E.oS(1):E.bS(1):E.dS(1);vA=a:min(a+E.bS(1)-1,E.dS(1));
    xT=x;
    xouT=xT;
    for b=E.oS(2):E.bS(2):E.dS(2);vB=b:min(b+E.bS(2)-1,E.dS(2));
        xS=dynInd(xT,vB,4);
        %FOURIER DOMAIN (SEGMENTS)
        xouS=xS;
        if isfield(E,'Fs') && ~isempty(E.Fs)
            for c=1:length(vA)
                xR=dynInd(xS,E.nEc(vA(c))+1:E.nEc(vA(c)+1),1);
                if isfield(EH,'We') && ~isempty(EH.We);xR=bsxfun(@times,xR,EH.We(E.nEc(vA(c))+1:E.nEc(vA(c)+1)));end
                xouR=[];
                if ~isempty(E.Fs{1}{vA(c)})
                    Nec=size(E.Fs{1}{vA(c)},1);
                    for d=1:inBlSz:Nec;vD=d:min(d+inBlSz-1,Nec);                     
                        if size(E.Fs,1)==2;xQ=bsxfun(@times,dynInd(xR,vD,1),conj(dynInd(E.Fs{2}{vA(c)},vD,1)));end  
                        xQ=aplGPU(dynInd(E.Fs{1}{vA(c)},vD,1)',xQ,1);
                        if d==1;xouR=xQ;else xouR=xouR+xQ;end
                    end
                end;xQ=[];
                if c==1 || isempty(xouS);xouS=xouR;elseif ~isempty(xouR);xouS=cat(5,xouS,xouR);end
            end;xR=[];xouR=[];
        end
        
        if isfield(E,'Es') && E.Es==1 && isfield(E,'Ef');xouS=aplGPU(E.Ef,xouS,E.pe);end%Fourier       
        %SENSE
        if isfield(EH,'Ub')%Sense unfolding (third dimension)             
            for n=3
                if isfield(EH,'Bb');xouS=repmat(xouS,[1 1 EH.Ub{n}.NX/EH.Ub{n}.NY]);elseif ~isempty(EH.Ub{n});xouS=ifold(xouS,n,EH.Ub{n}.NX,EH.Ub{n}.NY,EH.UAb{n});end
            end   
        end
        
        %MB/NYQUIST
        if isfield(EH,'Gh');xouS=bsxfun(@times,xouS,EH.Gh.Ab);end%Nyquist ghosting decoding
        if isfield(EH,'Bb');xouS=bsxfun(@times,xouS,dynInd(EH.Bb,E.cc,6));end%MB decoding
        if isfield(EH,'Eb');xouS=aplGPU(EH.Eb,xouS,E.pe);end%Space
        
        if isfield(EH,'Ub')%Sense unfolding (first two dimensions)
            for n=2:-1:1
                if ~isempty(EH.Ub{n});xouS=ifold(xouS,n,EH.Ub{n}.NX,EH.Ub{n}.NY,EH.UAb{n});end
            end
        end
        %COIL PROFILES
        if isfield(E,'Sf')
            Saux=dynInd(E.Sf,vB,4);
            if isa(xouS,'gpuArray');Saux=gpuArray(Saux);end
            xouS=bsxfun(@times,xouS,conj(Saux));
            if isfield(E,'Zf') && ~isempty(E.Zf);xouS=bsxfun(@times,xouS,dynInd(E.Zf,vB,4));end
            if ~isfield(E,'Zf') || ~isempty(E.Zf);xouS=sum(xouS,4);end
        end%Sensitivities
        if b==E.oS(2) || isempty(xouT);xouT=xouS;elseif isfield(E,'Sf') && ~isempty(xouS);xouT=xouT+xouS;elseif ~isempty(xouS);xouT=cat(4,xouT,xouS);end
    end;xS=[];xouS=[];
    %RIGID TRANSFORM (MOTION STATES)
    if isfield(E,'Tr') && ~isempty(E.Tr)    
        if any(E.Tr(:)~=0) && any(vA<=E.NMs)
            if isfield(EH,'Tb');Tb=extractFactorsSincRigidTransform(EH.Tb,vA(vA<=E.NMs),5);
            else Tb=precomputeFactorsSincRigidTransform(E.kG,E.rkG,dynInd(E.Tr,vA(vA<=E.NMs),5),0,0,1,1);
            end
            xouT=dynInd(xouT,vA<=E.NMs,5,sincRigidTransform(dynInd(xouT,vA<=E.NMs,5),Tb,0,E.Fof,E.Fob,0));
        end
        if ~isempty(xouT);xouT=sum(xouT,5);end
    end    
    if a==E.oS(1) || isempty(xou);xou=xouT;elseif isfield(E,'Tr') && ~isempty(xouT);xou=xou+xouT;elseif ~isempty(xouT);xou=cat(5,xou,xouT);end    
end;xT=[];xouT=[];
if isfield(EH,'Fb')%Filtering
    x=xou;   
    if isfield(EH.Fb,'Fd')
        for n=1:EH.Fb.ND
            padEl=zeros(1,EH.Fb.ND+1);padEl(n)=1;
            xs=diff(padarray(dynInd(x,n,EH.Fb.ND+1),padEl,0,'pre'),1,n)/EH.Fb.w(n);
            if n==1;xou=xs;else xou=xou+xs;end
        end            
    else
        ND=numDims(x);
        cont=1;
        if length(EH.Fb)>1;extDim=1;else extDim=0;end
        for n=1:length(EH.Fb)        
            mirr=zeros(1,ND-extDim);mirr(n)=EH.ma(n);       
            if EH.ma(n)
                xs=cat(n,dynInd(x,cont,ND),dynInd(x,cont+1,ND));cont=cont+2;
            else
                xs=dynInd(x,cont,ND+1-extDim);cont=cont+1;
            end
            if isa(EH.Fb{n},'gpuArray');xs=gpuArray(xs);end
            xs=filtering(xs,EH.Fb{n},EH.mi(n));
            xs=mirroring(xs,mirr,0);
            if n==1;xou=xs;else xou=xou+xs;end
        end
    end
end

if isfield(EH,'NXRec');xou=resampling(xou,EH.NXRec);end

if isfield(EH,'Mc');xou=bsxfun(@times,xou,EH.Mc);end%Mask for sensitivity computation
if isfield(EH,'Xb');xou=bsxfun(@times,xou,EH.Xb);end%Body coil or mask
if isfield(EH,'Gb');xou=filtering(xou,EH.Gb,EH.mi);end%Filtering

if isfield(EH,'Ab');xou=EH.Ab*xou;end%Explicit matrix


