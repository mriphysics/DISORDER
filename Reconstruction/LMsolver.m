function E=LMsolver(y,E,x,C,deb)

%LMSOLVER   Performs a Levenger-Marquardt iteration on a series of 
%   parameters of the encoding matrix to minimize a least squares 
%   backprojection of the reconstruction to the measured data
%   E=LMsolver(Y,E,X,{C},{DEB})
%   * Y is the measured data
%   * E is the encoding structure
%   * X is the reconstructed data
%   * {C} is a constrain structure
%   * {DEB} indicates whether to print information about convergence
%   ** E is the modified encoding structure
%

if nargin<4;C=[];end
if nargin<5 || isempty(deb);deb=2;end

if isfield(E,'Tr')              
    %GENERAL PARAMETERS AND ARRAYS
    gpu=isa(x,'gpuArray');
    T=E.Tr;E=rmfield(E,'Tr');
    NT=size(T);
    NX=size(x);    
    dimG=max(6,length(NT));dimS=find(NT(1:dimG-1)~=1);if isempty(dimS);dimS=dimG-1;end%Dimensions of motion states and parameters of the transform 
    ndS=NT(dimS);ndG=NT(dimG);NTS=length(dimS);
    ha=horzcat(repmat(1:ndG,[2 1]),nchoosek(1:ndG,2)');%Combinations of derivatives with repetitions to approximate the Hessian
    ndH=size(ha,2);   
    flagw=zeros(NT(1:dimS(end)));
    En=single(zeros([ndS 1]));EnPrev=En;EnPrevF=EnPrev;EnF=En;   
    dH=single(zeros([ndS ndH]));dG=single(zeros([ndS ndG]));dGEff=dG;        
    multA=1.2;multB=2;%Factors to divide/multiply the weight that regularizes the Hessian matrix when E(end)<E(end-1)    
    NElY=numel(y);
    if isfield(E,'nF');E.Sf=dynInd(E.Sf,E.nF,3);y=dynInd(y,E.nF,3);end%Extract ROI in the readout direction
    if isfield(E,'Ps') && ~isempty(E.Ps);E.Sf=bsxfun(@times,E.Sf,E.Ps);y=bsxfun(@times,y,E.Ps);end%Preconditioner of the coils, deprecated
    
    %BLOCK SIZES AND CONVERGENCE VALUES FOR THE TRANSFORM
    ET.bS=E.bS;ET.dS=E.dS;ET.oS=E.oS;ET.kG=E.kG;ET.rkG=E.rkG;   
    if isfield(E,'Bm');ET.Bm=E.Bm;ET.bS(1)=1;end
    if isfield(E,'Fms');ET.Sf=E.Sf;ET.Fms=E.Fms;ET.bS(1)=NX(3);E.dS(1)=1;end
    if isfield(E,'ZSl');ET.ZSl=E.ZSl;E.ZSl=-E.ZSl;end
    ET.cTFull=E.cT;       

    if NTS==1;NB=1;else NB=NT(dimS(1));end
    for b=1:NB
        if NB==1;ET.cT=ET.cTFull;else ET.cT=dynInd(ET.cTFull,b,dimS(1));end
        if isfield(E,'Fms');E.Fms=dynInd(ET.Fms,b,5);end

        %COMPUTE THE JACOBIAN AND STORE VALUES
        ind2EstOr=find(~ET.cT(:));NEst=length(ind2EstOr);
        if isfield(E,'Fs');ry=zeros([size(y,1) 1],'like',real(y));else ry=zeros([NT(6) 1],'like',real(y));end

        for a=1:ET.bS(1):NEst;vA=a:min(a+ET.bS(1)-1,NEst);vA=ind2EstOr(vA);
            if isfield(E,'Fms')
               [~,E.kG,E.rkG]=generateTransformGrids(NX,gpu,NX,ceil((NX+1)/2),1,[],abs(E.ZSl),vA);
               E.Sf=dynInd(ET.Sf,vA,3);            
            end                               
            if isfield(E,'Fs');indY=[];for c=1:length(vA);indY=horzcat(indY,E.nSt(vA(c))+1:E.nSt(vA(c)+1));end;end
            if NB==1;vT=vA;else vT={b,vA};end
            Ti=dynInd(T,vT,dimS);
            [Tf,Tg]=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Ti,1,1,1,1);%Transform factors 

            xT=x;
            if isfield(E,'Dc');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,E.Dc.d,6),E.Dc.D));end
            if isfield(E,'ZSl');xT=extractSlabs(xT,abs(E.ZSl),1,1);end
            if isfield(E,'Fms');xT=dynInd(xT,vA,6);end                       
            [xT,xTG]=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);
            
            if isfield(E,'nF');xT=dynInd(xT,E.nF,3);end
            if isfield(E,'Fs')
                E.vA=vA;
                E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
                xT=encode(xT,E)-dynInd(y,indY,1);
            elseif isfield(E,'Bm')
                E.Bm=dynInd(ET.Bm,vA,6);
                xT=encode(xT,E)-bsxfun(@times,E.Bm,y);
            else
                xT=encode(xT,E)-dynInd(y,vT,[5 3]);
            end
            
            if isfield(E,'Pf') && ~isempty(E.Pf);xT=bsxfun(@times,xT,dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
            if isfield(E,'Fs')
                ry(indY)=normm(xT,[],3:4)/NElY;
                EnPrevF=dynInd(EnPrevF,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),ry(indY)),vA,1)));
            else
                ry(vA)=normm(xT)/NElY;
                EnPrevF=dynInd(EnPrevF,vT,1:NTS,gather(ry(vA)));
            end
            if isfield(E,'Fd') && ~isempty(E.Fd)%Filtering
                xT=fftGPU(xT,3);
                xT=bsxfun(@times,xT,dynInd(E.Fd,indY,1));
                EnPrev=dynInd(EnPrev,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),normm(xT,[],3:4)/NElY),vA,1)));
            else
                EnPrev=dynInd(EnPrev,vT,1:NTS,dynInd(EnPrevF,vT,1:NTS));
            end

            xTG=sincRigidTransformGradient(xTG,Tf,Tg,E.Fof,E.Fob);
            for g=1:ndG
                if isfield(E,'nF');xTG{g}=dynInd(xTG{g},E.nF,3);end
                xTG{g}=encode(xTG{g},E);
                if isfield(E,'Pf') && ~isempty(E.Pf);xTG{g}=bsxfun(@times,xTG{g},dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
                if isfield(E,'Fd') && ~isempty(E.Fd)
                    xTG{g}=fftGPU(xTG{g},3);
                    xTG{g}=bsxfun(@times,xTG{g},dynInd(E.Fd,indY,1));
                end
                if isfield(E,'Fs')
                    dG(vA,g)=gather(dynInd(accumarray(E.mSt(indY),multDimSum(real(xTG{g}.*conj(xT)),3:4)),vA,1));
                elseif isfield(E,'Bm')
                    dG(vA,g)=gather(multDimSum(real(xTG{g}.*conj(xT)),1:6));
                else
                    dG=dynInd(dG,[vT g],1:NTS+1,gather(permute(multDimSum(real(xTG{g}.*conj(xT)),[1:2 4]),[5 3 1 2 4])));
                end
            end
            for h=1:ndH
                if isfield(E,'Fs')
                    dH(vA,h)=gather(dynInd(accumarray(E.mSt(indY),multDimSum(real(xTG{ha(1,h)}.*conj(xTG{ha(2,h)})),3:4)),vA,1));
                elseif isfield(E,'Bm')
                    dH(vA,h)=gather(multDimSum(real(xTG{ha(1,h)}.*conj(xTG{ha(2,h)})),1:6));
                else
                    dH=dynInd(dH,[vT h],1:NTS+1,gather(permute(multDimSum(real(xTG{ha(1,h)}.*conj(xTG{ha(2,h)})),[1:2 4]),[5 3 1 2 4])));
                end
            end
        end;xTG=[];
    end

    %UPDATE
    En=EnPrev;EnF=EnPrevF;
    MH=single(eye(ndG));
    fina=0;
    while ~fina   
        %BUILD HESSIAN MATRIX AND POTENTIAL UPDATE
        ry(:)=0;
        for a=find(~ET.cTFull(:) & E.w(:)<1e10)'            
            aa=ind2subV(NT(dimS),a);
            for h=1:ndH
                if ha(1,h)==ha(2,h);MH(ha(1,h),ha(2,h))=(1+E.w(a))*dynInd(dH,[aa h],1:NTS+1);else MH(ha(1,h),ha(2,h))=dynInd(dH,[aa h],1:NTS+1);MH(ha(2,h),ha(1,h))=dynInd(dH,[aa h],1:NTS+1);end
            end      
            GJ=double(dynInd(dG,aa,1:NTS));
            GJ=single(double(MH)\GJ(:));
            GJ=resPop(GJ,1,[],NTS+1);
            dGEff=dynInd(dGEff,aa,1:NTS,-(E.winit/E.w(a))*GJ);
        end  
        Tupr=shiftdim(dGEff,-(dimS(1)-1));
        Tup=restrictTransform(T+Tupr);%Rotation parameters between -pi and pi

        %CHECK ENERGY REDUCTION
        E.cTFull=(ET.cTFull | flagw);
        for b=1:NB
            if NB==1;E.cT=E.cTFull;else E.cT=dynInd(E.cTFull,b,dimS(1));end
            if isfield(E,'Fms');E.Fms=dynInd(ET.Fms,b,5);end
                                        
            ind2Est=find(~E.cT(:));NEst=length(ind2Est);
            for a=1:ET.bS(1):NEst;vA=a:min(a+ET.bS(1)-1,NEst);vA=ind2Est(vA);
                if isfield(E,'Fms')
                    [~,E.kG,E.rkG]=generateTransformGrids(NX,gpu,NX,ceil((NX+1)/2),1,[],abs(E.ZSl),vA);
                    E.Sf=dynInd(ET.Sf,vA,3);                   
                end                
                if isfield(E,'Fs');indY=[];for c=1:length(vA);indY=horzcat(indY,E.nSt(vA(c))+1:E.nSt(vA(c)+1));end;end
                if NB==1;vT=vA;else vT={b,vA};end
                Ti=dynInd(Tup,vT,dimS);
                Tf=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Ti,1,0,1,1);

                xT=x;               
                if isfield(E,'Dc');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,E.Dc.d,6),E.Dc.D));end
                if isfield(E,'ZSl');xT=extractSlabs(xT,abs(E.ZSl),1,1);end
                if isfield(E,'Fms');xT=dynInd(xT,vA,6);end                
                xT=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);   
                if isfield(E,'nF');xT=dynInd(xT,E.nF,3);end            
                if isfield(E,'Fs')
                    E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
                    E.vA=vA;
                    xT=encode(xT,E)-dynInd(y,indY,1);
                elseif isfield(E,'Bm')
                    E.Bm=dynInd(ET.Bm,vA,6);
                    xT=encode(xT,E)-bsxfun(@times,E.Bm,y);
                else
                    xT=encode(xT,E)-dynInd(y,vT,[5 3]);
                end
                if isfield(E,'Pf') && ~isempty(E.Pf);xT=bsxfun(@times,xT,dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
                if isfield(E,'Fs')
                    ry(indY)=normm(xT,[],3:4)/NElY;
                    EnF=dynInd(EnF,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),ry(indY)),vA,1)));
                else
                    ry(vA)=normm(xT)/NElY;
                    EnF=dynInd(EnF,vT,1:NTS,gather(ry(vA)));                                 
                end
                if isfield(E,'Fd') && ~isempty(E.Fd)
                    xT=fftGPU(xT,3);
                    xT=bsxfun(@times,xT,dynInd(E.Fd,indY,1));
                    En=dynInd(En,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),normm(xT,[],3:4)/NElY),vA,1)));
                else
                    En=dynInd(En,vT,1:NTS,dynInd(EnF,vT,1:NTS));
                end   
            end
        end
        En(E.w(:)>1e10)=EnPrev(E.w(:)>1e10);         
        flagw(En(:)<=EnPrev(:))=1;

        if deb>=2;fprintf('Energy before: %0.6g / Energy after: %0.6g\n',sum(EnPrevF(:)),sum(EnF(:)));end  
        if any(~flagw(:)) 
            E.w(En(:)>EnPrev(:) & ~ET.cTFull(:))=E.w(En(:)>EnPrev(:) & ~ET.cTFull(:))*multB;
        else
            E.w(~ET.cTFull(:))=E.w(~ET.cTFull(:))/multA;
            E.w(E.w<1e-8)=multA*E.w(E.w<1e-8);%To avoid numeric instabilities 
            T=Tup;
            fina=1;
            traDiff=abs(dynInd(Tupr,1:3,dimG));
            rotDiff=abs(convertRotation(dynInd(Tupr,4:6,dimG),'rad','deg'));  
            traDiffMax=multDimMax(traDiff,dimS);rotDiffMax=multDimMax(rotDiff,dimS);    
            ET.cTFull(max(traDiff,[],dimG)<E.tL(1) & max(rotDiff,[],dimG)<E.tL(2))=1;
            if deb>=2
                fprintf('Maximum change in translation (vox): %s/ rotation (deg): %s\n',sprintf('%0.3f ',traDiffMax(:)),sprintf('%0.3f ',rotDiffMax(:)));
                fprintf('Not converged motion states: %d of %d\n',prod(ndS)-sum(single(ET.cTFull(:))),prod(ndS));
            end              
        end
    end  
        
    T=constrain(T,C);
    E.Tr=T;E.cT=ET.cTFull;E.bS=ET.bS;E.dS=ET.dS;E.oS=ET.oS;E.kG=ET.kG;E.rkG=ET.rkG;
    if isfield(E,'Fms');E.Sf=ET.Sf;E.Fms=ET.Fms;end
    if isfield(E,'ZSl');E.ZSl=ET.ZSl;end
    if isfield(E,'Bm');E.Bm=ET.Bm;end
    if isfield(E,'vA');E=rmfield(E,'vA');end
    if isfield(E,'Fd');E=rmfield(E,'Fd');end   
    Enmin=min(EnPrev,En);
    if ~isfield(E,'En');E.En=Enmin;else E.En(ind2EstOr)=Enmin(ind2EstOr);end%To keep a record of the achieved energy
end
