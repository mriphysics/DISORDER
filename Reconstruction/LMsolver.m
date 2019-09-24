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

gpu=isa(y,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

if isfield(E,'Tr')              
    %GENERAL PARAMETERS AND ARRAYS
    T=E.Tr;E=rmfield(E,'Tr');
    dimS=5;dimG=6;%Dimensions of motion states and parameters of the transform
    NT=size(T);ndS=NT(dimS);ndG=NT(dimG);        
    ha=horzcat(repmat(1:ndG,[2 1]),nchoosek(1:ndG,2)');%Combinations of derivatives with repetitions to approximate the Hessian
    ndH=size(ha,2);   
    flagw=zeros(NT(1:dimS));
    En=single(zeros([ndS 1]));EnPrev=En;EnPrevF=EnPrev;EnF=En;   
    dH=single(zeros([ndS ndH]));dG=single(zeros([ndS ndG]));dGEff=dG;        
    multA=1.2;multB=2;%Factors to divide/multiply the weight that regularizes the Hessian matrix when E(end)<E(end-1)    
    NElY=numel(y);
    if isfield(E,'nF');E.Sf=dynInd(E.Sf,E.nF,3);y=dynInd(y,E.nF,3);end%Extract ROI in the readout direction    
    if isfield(E,'Ps') && ~isempty(E.Ps);E.Sf=bsxfun(@times,E.Sf,E.Ps);y=bsxfun(@times,y,E.Ps);end%Preconditioner of the coils, deprecated
    
    %BLOCK SIZES AND CONVERGENCE VALUES FOR THE TRANSFORM
    ET.bS=E.bS;ET.dS=E.dS;ET.oS=E.oS;
    ET.cT=E.cT;            
    
    %COMPUTE THE JACOBIAN AND STORE VALUES
    ind2EstOr=find(~ET.cT(:));NEst=length(ind2EstOr);
    ry=zeros([size(y,1) 1],'like',real(y));          
    for a=1:ET.bS(1):NEst;vA=a:min(a+ET.bS(1)-1,NEst);vA=ind2EstOr(vA);
        indY=[];for c=1:length(vA);indY=horzcat(indY,E.nSt(vA(c))+1:E.nSt(vA(c)+1));end
        [Tf,Tg]=precomputeFactorsSincRigidTransform(E.kG,E.rkG,dynInd(T,vA,dimS),1,1,1,1);%Transform factors        
        
        [xT,xTG]=sincRigidTransform(x,Tf,1,E.Fof,E.Fob);
        if isfield(E,'nF');xT=dynInd(xT,E.nF,3);end
        E.vA=vA;
        E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
        xT=encode(xT,E)-dynInd(y,indY,1);       
        if isfield(E,'Pf') && ~isempty(E.Pf);xT=bsxfun(@times,xT,dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
        ry(indY)=normm(xT,[],3:4)/NElY;
        EnPrevF(vA)=gather(dynInd(accumarray(E.mSt(indY),ry(indY)),vA,1));
        if isfield(E,'Fd') && ~isempty(E.Fd)%Filtering
            xT=fftGPU(xT,3,gpuF);
            xT=bsxfun(@times,xT,dynInd(E.Fd,indY,1));
            EnPrev(vA)=gather(dynInd(accumarray(E.mSt(indY),normm(xT,[],3:4)/NElY),vA,1));
        else
            EnPrev(vA)=EnPrevF(vA);
        end
  
        xTG=sincRigidTransformGradient(xTG,Tf,Tg,E.Fof,E.Fob);
        for g=1:ndG
            if isfield(E,'nF');xTG{g}=dynInd(xTG{g},E.nF,3);end
            xTG{g}=encode(xTG{g},E);
            if isfield(E,'Pf') && ~isempty(E.Pf);xTG{g}=bsxfun(@times,xTG{g},dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
            if isfield(E,'Fd') && ~isempty(E.Fd)
                xTG{g}=fftGPU(xTG{g},3,gpuF);
                xTG{g}=bsxfun(@times,xTG{g},dynInd(E.Fd,indY,1));
            end
            dG(vA,g)=gather(dynInd(accumarray(E.mSt(indY),multDimSum(real(xTG{g}.*conj(xT)),3:4)),vA,1));
        end
        for h=1:ndH;dH(vA,h)=gather(dynInd(accumarray(E.mSt(indY),multDimSum(real(xTG{ha(1,h)}.*conj(xTG{ha(2,h)})),3:4)),vA,1));end
    end;xTG=[];
    
    %UPDATE
    En=EnPrev;EnF=EnPrevF;
    MH=single(eye(ndG));        
    fina=0;
    while ~fina   
        %BUILD HESSIAN MATRIX AND POTENTIAL UPDATE
        ry(:)=0;
        for a=find(~ET.cT & E.w<1e10)'           
            for h=1:ndH
                if ha(1,h)==ha(2,h);MH(ha(1,h),ha(2,h))=(1+E.w(a))*dH(a,h);else MH(ha(1,h),ha(2,h))=dH(a,h);MH(ha(2,h),ha(1,h))=dH(a,h);end
            end      
            dGEff(a,:)=-(E.winit/E.w(a))*single(double(MH)\double(dG(a,:)'))';         
        end  
        Tupr=shiftdim(dGEff,-(dimS-1));
        Tup=restrictTransform(T+Tupr);%Rotation parameters between -pi and pi
        
        %CHECK ENERGY REDUCTION
        E.cT=(ET.cT | flagw);
        ind2Est=find(~E.cT(:));NEst=length(ind2Est);        
        for a=1:ET.bS(1):NEst;vA=a:min(a+ET.bS(1)-1,NEst);vA=ind2Est(vA);        
            indY=[];for c=1:length(vA);indY=horzcat(indY,E.nSt(vA(c))+1:E.nSt(vA(c)+1));end
            Tf=precomputeFactorsSincRigidTransform(E.kG,E.rkG,dynInd(Tup,vA,dimS),1,0,1,1);
            xT=sincRigidTransform(x,Tf,1,E.Fof,E.Fob);   
            if isfield(E,'nF');xT=dynInd(xT,E.nF,3);end            
            E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
            E.vA=vA;
            xT=encode(xT,E)-dynInd(y,indY,1);
            if isfield(E,'Pf') && ~isempty(E.Pf);xT=bsxfun(@times,xT,dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
            ry(indY)=normm(xT,[],3:4)/NElY;
            EnF(vA)=gather(dynInd(accumarray(E.mSt(indY),ry(indY)),vA,1));
            if isfield(E,'Fd') && ~isempty(E.Fd)
                xT=fftGPU(xT,3,gpuF);
                xT=bsxfun(@times,xT,dynInd(E.Fd,indY,1));
                En(vA)=gather(dynInd(accumarray(E.mSt(indY),normm(xT,[],3:4)/NElY),vA,1));
            else
                En(vA)=EnF(vA);
            end   
        end
        En(E.w(:)>1e10)=EnPrev(E.w(:)>1e10);            
        flagw(En<=EnPrev)=1;
        
        if deb>=2;fprintf('Energy before: %0.6g / Energy after: %0.6g\n',sum(EnPrevF),sum(EnF));end  
        if any(~flagw) 
            E.w(shiftdim(En>EnPrev,-(dimS-1)) & ~ET.cT)=E.w(shiftdim(En>EnPrev,-(dimS-1)) & ~ET.cT)*multB;
        else
            E.w(~ET.cT)=E.w(~ET.cT)/multA;
            E.w(E.w<1e-8)=multA*E.w(E.w<1e-8);%To avoid numeric instabilities 
            T=Tup;
            fina=1;
            traDiff=abs(dynInd(Tupr,1:3,dimG));
            rotDiff=abs(convertRotation(dynInd(Tupr,4:6,dimG),'rad','deg'));  
            traDiffMax=max(traDiff,[],dimS);rotDiffMax=max(rotDiff,[],dimS);    
            ET.cT(max(traDiff,[],dimG)<E.tL(1) & max(rotDiff,[],dimG)<E.tL(2))=1;
            if deb>=2
                fprintf('Maximum change in translation (vox): %s/ rotation (deg): %s\n',sprintf('%0.3f ',traDiffMax(:)),sprintf('%0.3f ',rotDiffMax(:)));
                fprintf('Not converged motion states: %d of %d\n',ndS-sum(single(ET.cT)),ndS);
            end              
        end       
    end
    T=constrain(T,C);
    E.Tr=T;E.cT=ET.cT;E.bS=ET.bS;E.dS=ET.dS;E.oS=ET.oS;
    if isfield(E,'vA');E=rmfield(E,'vA');end
    if isfield(E,'Fd');E=rmfield(E,'Fd');end   
    Enmin=min(EnPrev,En);
    if ~isfield(E,'En');E.En=Enmin;else E.En(ind2EstOr)=Enmin(ind2EstOr);end%To keep a record of the achieved energy
end
