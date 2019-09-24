function xou=encodedecode(x,E,EH)

%ENCODEDECODE   Applies a given encoding+decoding to the data
%   XOU=ENCODEDECODE(X,E)
%   * X is the data before encoding+decoding
%   * E is the encoding structure
%   * EH is the decoding structure
%   ** XOU is the data after encoding+decoding
%

if isempty(E);xou=x;return;end
inBlSz=1e3;

if ~isfield(E,'bS')%Block sizes by default
    E.bS=[size(x,5) size(x,4)];E.dS=[size(x,5) size(x,4)];E.oS=[1 1];
    if isfield(E,'Sf');E.bS(2)=size(E.Sf,4);E.dS(2)=size(E.Sf,4);end
    if isfield(E,'Tr');E.bS(1)=E.NSe;E.dS(1)=E.NSe;end
end
if isfield(E,'Zf') && isempty(E.Zf) && isfield(E,'Zf');E.oS(2)=1;E.bS(2)=size(E.Sf,4);E.dS(2)=size(E.Sf,4);end

if isfield(E,'Af');x=E.Af*x;end%Explicit matrix

if isfield(E,'Ti');xTi=E.Ti.la*x;end%This does not have counterpart in decoding
if isfield(E,'Gf');x=filtering(x,E.Gf,E.mi);end%Filtering
if isfield(E,'Xf');x=bsxfun(@times,x,E.Xf);end%Body coil or mask
if isfield(E,'Ti');x=x+xTi;end%This does not have counterpart in decoding
xou=x;
if isfield(E,'Ff')%Filtering
    xou=[];
    if isfield(E.Ff,'Fd')
        for n=1:E.Ff.ND
            padEl=zeros(1,E.Ff.ND+1);padEl(n)=1;
            xs=padarray(diff(x,1,n)/E.Ff.w(n),padEl,0,'post');
            if n==1;xou=xs;else xou=cat(E.Ff.ND+1,xou,xs);end  
        end
    else
        ND=numDims(x);    
        N=size(x);
        for n=1:length(E.Ff)
            mirr=zeros(1,ND);mirr(n)=E.ma(n);
            xs=mirroring(x,mirr,1);
            xs=filtering(xs,E.Ff{n},E.mi(n));
            xs=resPop(xs,n,[N(n) 1+E.ma(n)],[n ND+1]);
            if E.Sm;xs=gather(xs);end
            if n==1;xou=xs;else xou=cat(ND+1,xou,xs);end
        end
    end
    x=xou;
end

if isfield(E,'NXAcq');x=resampling(x,E.NXAcq);end

for a=E.oS(1):E.bS(1):E.dS(1);vA=a:min(a+E.bS(1)-1,E.dS(1));
    if isfield(E,'vA');vA=vA(ismember(vA,E.vA));end
    xT=x;
    %RIGID TRANSFORM (MOTION STATES)
    if isfield(E,'Tr') && ~isempty(E.Tr)
        if any(vA<=E.NMs)
            if any(E.Tr(:)~=0)
                if isfield(E,'Tf');Tf=extractFactorsSincRigidTransform(E.Tf,vA(vA<=E.NMs),5);
                else Tf=precomputeFactorsSincRigidTransform(E.kG,E.rkG,dynInd(E.Tr,vA(vA<=E.NMs),5),1,0,1,1);
                end
                xT=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);
            else
                xT=repmat(xT,[ones(1,4) length(vA(vA<=E.NMs))]);
            end
        end
        if any(vA>E.NMs)
            assert(~isempty(xT),'Empty transformed array');
            xT=cat(5,xT,repmat(x,[ones(1,4) sum(vA>E.NMs)]));
        end
    end      
    xouT=xT;
    %COIL PROFILES
    for b=E.oS(2):E.bS(2):E.dS(2);vB=b:min(b+E.bS(2)-1,E.dS(2));
        xS=xT;
        if isfield(E,'Sf')
            Saux=dynInd(E.Sf,vB,4);
            if isa(xS,'gpuArray');Saux=gpuArray(Saux);end
            xS=bsxfun(@times,xS,Saux);
            xS=sum(xS,6);
        end%Sensitivities
        %SENSE
        if isfield(E,'Uf')%Sense folding (first two dimensions)
            for n=1:2
                if ~isempty(E.Uf{n});xS=fold(xS,n,E.Uf{n}.NX,E.Uf{n}.NY,E.UAf{n});end
            end
        end
        
        %MB/NYQUIST
        if isfield(E,'Ef');xS=aplGPU(E.Ef,xS,2);end%Fourier
        if isfield(E,'Bf');xS=bsxfun(@times,xS,dynInd(E.Bf,E.cc,6));end%MB encoding
        if isfield(E,'Gh');xS=bsxfun(@times,xS,E.Gh.Af);end%Nyquist ghosting encoding
        
        if isfield(E,'Uf')%Sense folding (third dimension)    
            for n=3
                if isfield(E,'Bf');xS=resSub(sum(resSub(xS,3,[E.Uf{n}.NY E.Uf{n}.NX/E.Uf{n}.NY]),4),3:4);elseif ~isempty(E.Uf{n});xS=fold(xS,n,E.Uf{n}.NX,E.Uf{n}.NY,E.UAf{n});end
            end
        end
        
        %FOURIER DOMAIN (SEGMENTS)
        xouS=xS;
        if isfield(E,'Fs') && ~isempty(E.Fs)
            for c=1:length(vA)
                xR=dynInd(xS,c,5);
                if ~isempty(E.Fs{1}{vA(c)})
                    xouR=xR;xouR(:)=0;
                    Nec=size(E.Fs{1}{vA(c)},1);
                    if isfield(EH,'We') && ~isempty(EH.We);vE=E.nEc(vA(c))+1:E.nEc(vA(c)+1);end
                    for d=1:inBlSz:Nec;vD=d:min(d+inBlSz-1,Nec);
                        xQ=xR;
                        F1=dynInd(E.Fs{1}{vA(c)},vD,1);
                        xQ=aplGPU(F1,xQ,1);
                        if size(E.Fs,1)==2
                            F2=dynInd(E.Fs{2}{vA(c)},vD,1);
                            xQ=sum(bsxfun(@times,xQ,F2),2);
                            if isfield(EH,'We') && ~isempty(EH.We);xQ=bsxfun(@times,xQ,EH.We(vE(vD)));end
                            xQ=bsxfun(@times,xQ,conj(F2));
                        end 
                        xouR=xouR+aplGPU(F1',xQ,1);
                    end;xQ=[]; 
                else
                    xouR=[];
                end
                if c==1 || isempty(xouS);xouS=xouR;elseif ~isempty(xouR);xouS=cat(5,xouS,xouR);end
            end;xR=[];xouR=[];
        end
        
        %SENSE
        if isfield(EH,'Ub')%Sense unfolding (third dimension)             
            for n=3        
                if isfield(EH,'Bb');xouS=repmat(xouS,[1 1 EH.Ub{n}.NX/EH.Ub{n}.NY]);elseif ~isempty(EH.Ub{n});xouS=ifold(xouS,n,EH.Ub{n}.NX,EH.Ub{n}.NY,EH.UAb{n});end      
            end   
        end
        
        %MB/NYQUIST
        if isfield(EH,'Gh');xouS=bsxfun(@times,xouS,EH.Gh.Ab);end%Nyquist ghosting decoding
        if isfield(EH,'Bb');xouS=bsxfun(@times,xouS,dynInd(EH.Bb,E.cc,6));end%MB decoding   
        if isfield(EH,'Eb');xouS=aplGPU(EH.Eb,xouS,2);end%Fourier
        
        if isfield(EH,'Ub')%Sense unfolding (first two dimensions)
            for n=2:-1:1
                if ~isempty(EH.Ub{n});xouS=ifold(xouS,n,EH.Ub{n}.NX,EH.Ub{n}.NY,EH.UAb{n});end
            end
        end
        %COIL PROFILES
        if isfield(E,'Sf')        
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

if isfield(EH,'NXRec');xou=resampling(xou,EH.NXRec);end

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

if isfield(EH,'Mc');xou=bsxfun(@times,xou,EH.Mc);end%Mask for sensitivity computation
if isfield(EH,'Xb');xou=bsxfun(@times,xou,EH.Xb);end%Body coil or mask
if isfield(EH,'Gb');xou=filtering(xou,EH.Gb,EH.mi);end%Filtering

if isfield(EH,'Ab');xou=EH.Ab*xou;end%Explicit matrix
