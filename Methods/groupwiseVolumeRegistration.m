function [x,T,y]=groupwiseVolumeRegistration(x,W,T,meanT,lev,fracOrd,parComp)

%GROUPWISEVOLUMEREGISTRATION   Per-volume groupwise rigid registration of dynamic studies
%   [X,T]=GROUPWISEVOLUMEREGISTRATION(X,{W},{T},{MEANT},{LEV},{FRACORD},{PARCOMP})
%   * X is the data to be registered
%   * {W} is a mask for the ROI to compute the metric
%   * {T} is the initial motion parameters (defaults to 0)
%   * {MEANT} is a flag to project the solution to the space where the sum 
%   of the transform parameters through the volumes is 0, defaults to 1
%   * {LEV} are the resolution levels for correction
%   * {FRACORD} is the fractional finite difference metric order
%   * {PARCOMP} serves to compute only a given subset of parameters
%   ** X is the registered data
%   ** T are the estimated motion parameters
%   ** Y is the average and std in time
%

NX=size(x);NX(end+1:8)=1;ND=max(numDims(x),3);

if nargin<2 || isempty(W);W=ones(NX(1:min(ND-1,3)),'like',x);end
if nargin<3 || isempty(T);T=single(zeros([ones(1,ND-1) NX(ND) 6]));end
NXV=NX(1:3);
if nargin<4 || isempty(meanT);meanT=1;end
if nargin<5 || isempty(lev);lev=[2 1];end
if nargin<6 || isempty(fracOrd);fracOrd=0;end
if nargin<7 || isempty(parComp);parComp=ones(1,6);end

gpu=isa(x,'gpuArray');
NT=size(T);ndT=ndims(T);
perm=1:ndT;perm([2 ndT])=[ndT 2];
parComp=permute(parComp,perm);

a=[1 2 3 1 1 2 1 2 3 1 2 3 1 2 3 4 4 5 4 5 6;
   1 2 3 2 3 3 4 4 4 5 5 5 6 6 6 5 6 6 4 5 6];
NHe=size(a,2);
dHe=single(zeros([NHe NT(ndT-1)]));
dH=single(zeros(NT([ndT ndT-1])));dHEff=dH;
Eprev=single(zeros(NT(1:ndT-1)));E=single(zeros(NT(1:ndT-1)));
convT=single(false(NT(1:ndT-1)));
if gpu;[Eprev,E,T,convT]=parUnaFun({Eprev,E,T,convT},@gpuArray);end
multA=1.2;multB=2;%Factors to divide/multiply the weight that regularizes the Hessian matrix when E(end)<E(end-1)

BlSz=10;
for l=1:length(lev) 
    fprintf('Resolution level: %d\n',lev(l));
    NXres=round(NXV/lev(l));
    [~,kGrid,rkGrid]=generateTransformGrids(NXV,gpu,NXres,[],1);
    [FT,FTH]=buildStandardDFTM(NXres,0,gpu);

    fina=0;
    winic=1;
    w=winic*ones(NT(1:ndT-1));
    flagw=zeros(NT(1:ndT-1));

    y=single(zeros([NXres NX(4:ndT-2)]));
    if gpu;y=gpuArray(y);end
    perm=1:ndT;perm([1 2 ndT-1 ndT])=[ndT-1 ndT 1 2];
    
    xRes=resampling(x,NXres);WRes=abs(resampling(W,NXres));    
    H=buildFilter(2*NXres,'tukeyIso',ones(1,3),gpu,0.1,1);
    xRes=filtering(xRes,H,1); 
    
    if fracOrd~=0;GRes=buildFilter(NXres(1:2),'FractionalFiniteDiscreteIso',NX(1:2)./NXres(1:2),gpu,fracOrd);end%For fractional finite difference motion estimation

    cont=0;
    NcontTest=1;
    %Iterations
    while fina~=2
        dHe(:)=0;dH(:)=0;    
        if mod(cont,NcontTest)==0
            y(:)=0; 
            for s=1:BlSz(1):NT(ndT-1);vS=s:min(s+BlSz(1)-1,NT(ndT-1));%This is the ''reconstruction'' step
                etInv=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,vS,ndT-1),0,0,[],1);
                y=y+sum(sincRigidTransform(dynInd(xRes,vS,ndT-1),etInv,0,FT,FTH,0),ndT-1)/NT(ndT-1);
            end      
            convT(:)=0;
            cont=0;NcontTest=NcontTest+1;        
        end    
        cont=cont+1;

        for s=1:BlSz(1):NT(ndT-1);vS=s:min(s+BlSz(1)-1,NT(ndT-1));vS=vS(~convT(vS));
            if ~isempty(vS)
                [et,etg]=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,vS,ndT-1),1,1,[],1);
                [xT,xB]=sincRigidTransform(y,et,1,FT,FTH);
                WT=abs(sincRigidTransform(WRes,et,1,FT,FTH));
                xT=xT-dynInd(xRes,vS,ndT-1);
                if fracOrd~=0;xT=filtering(xT,GRes);end
                xT=bsxfun(@times,xT,WT);%IT SHOULD BE SQRT WT REALLY
                Eprev(vS)=multDimSum(real(xT.*conj(xT)),1:ndT-2);       
                G=sincRigidTransformGradient(xB,et,etg,FT,FTH);
                GH=cell(1,NT(ndT));
                for m=1:NT(ndT)
                    if fracOrd~=0;G{m}=filtering(G{m},GRes);end
                    G{m}=bsxfun(@times,G{m},WT);
                    GH{m}=conj(G{m});
                end
                for m=1:NHe   
                    GGE=real(G{a(1,m)}.*GH{a(2,m)});    
                    dHe(m,vS)=gather(permute(multDimSum(GGE,1:ndT-2),perm)); 
                end   
                for m=1:NT(ndT)       
                    G{m}=real(GH{m}.*xT);         
                    dH(m,vS)=gather(permute(multDimSum(G{m},1:ndT-2),perm));
                end
            end
        end
        
        %%%THERE'S A POTENTIAL ISSUE HERE IN THAT WE MAY KEEP ON UPDATING MOTION STATES THAT HAVE ALREADY CONVERGED. THIS MAY HAVE A SMALL EFFECT, THOUGH       
        MHe=single(eye(NT(ndT)));
        flagw(:)=0;    
        fina=0;
        while fina==0
            for s=1:NT(ndT-1)
                if ~convT(s)
                    for k=1:NHe
                        if a(1,k)==a(2,k);MHe(a(1,k),a(2,k))=(1+w(s))*dHe(k,s);
                        else MHe(a(1,k),a(2,k))=dHe(k,s);MHe(a(2,k),a(1,k))=dHe(k,s);
                        end              
                    end   
                    dHEff(:,s)=-winic*single(double(MHe)\double(dH(:,s)))/w(s);
                end
            end     
            dHEff(:,w>1e10)=0;
            permH=1:ndT;permH(ndT-1:ndT)=[2 1];permH(1:ndT-2)=3:ndT;
            Tupr=bsxfun(@times,parComp,permute(dHEff,permH));
            Tup=T+Tupr; 
            Tup=restrictTransform(Tup);                                             

            for s=1:BlSz(1):NT(ndT-1);vS=s:min(s+BlSz(1)-1,NT(ndT-1));vS=vS(~convT(vS) & flagw(vS)~=2);
                if ~isempty(vS)
                    et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(Tup,vS,ndT-1),1,0,[],1);                   
                    xT=sincRigidTransform(y,et,1,FT,FTH);
                    WT=abs(sincRigidTransform(WRes,et,1,FT,FTH));
                    xT=xT-dynInd(xRes,vS,ndT-1);
                    if fracOrd~=0;xT=filtering(xT,GRes);end
                    xT=bsxfun(@times,xT,WT);
                    E(vS)=multDimSum(real(xT.*conj(xT)),1:ndT-2);
                end
            end

            E(w>1e10)=Eprev(w>1e10);
            flagw(E<=Eprev)=2;               
            fprintf('Energy before: %0.6g / Energy after: %0.6g\n',sum(Eprev),sum(E));
            if any(flagw==1 | flagw==0)  
                w(E>Eprev & ~convT)=w(E>Eprev & ~convT)*multB;
            else   
                w(~convT)=w(~convT)/multA;
                w(w<1e-8)=multA*w(w<1e-8);%To avoid numeric instabilities 
                T=Tup;
                fina=2;
                traMax=abs(permute(dynInd(Tupr,1:3,ndT),perm));
                rotMax=convertRotation(abs(permute(dynInd(Tupr,4:6,ndT),perm)),'rad','deg');
                fprintf('Maximum change in translation (vox): ');fprintf('%0.3f ',max(traMax,[],1));
                fprintf('/ Maximum change in rotation (deg): ');fprintf('%0.3f ',max(rotMax,[],1));fprintf('\n');                       
                traLim=0.16;rotLim=0.08;                
                acce=lev(l);
                traLim=traLim*acce;rotLim=rotLim*acce;
                if max(traMax(:))>traLim || max(rotMax(:))>rotLim;fina=1;end                
                convT(max(traMax,[],2)<traLim & max(rotMax,[],2)<rotLim)=1;
                fprintf('Not converged motion states: %d of %d\n',NT(ndT-1)-sum(single(convT)),NT(ndT-1));
            end
        end
        if meanT;T=bsxfun(@minus,T,mean(T,ndT-1));end
    end
end

[~,kGrid,rkGrid]=generateTransformGrids(NXV,gpu,[],[],1);
[FT,FTH]=buildStandardDFTM(NXV,0,gpu);

for s=1:BlSz(1):NT(ndT-1);vS=s:min(s+BlSz(1)-1,NT(ndT-1));
    etInv=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,vS,ndT-1),0,0,[],1);
    x=dynInd(x,vS,ndT-1,sincRigidTransform(dynInd(x,vS,ndT-1),etInv,0,FT,FTH,0));
end
y=mean(x,ndT-1);
y=cat(ndT-1,y,std(x,0,ndT-1));
