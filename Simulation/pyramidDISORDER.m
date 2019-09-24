for s=1:size(NS,2)%Shots      
    for v=1:size(theta,2)%Motion levels
        for e=1:size(recTyp,1)%Reconstruction types
            if debug>0;fprintf('Solving for acquisition%s, motion%s, type %s\n',sprintf(' %d',NS(:,s)'),sprintf(' %.0f',theta(:,v)'),cell2mat(recTyp(e,:)));end
            
            %INITIALIZATION
            T=zeros(size(TGT{s}{v}),'like',TGT{s}{v});              
            Lk=str2double(recTyp(e,2));%Number of spatial multirresolution levels
            Lt=str2double(recTyp(e,3));%Number of temporal multirresolution levels
            nItPlan=round(cumsum(nExtern*ones(1,Lk+Lt-1)/(Lk+Lt-1)));%Iteration plan
            contKT=0;%Level counter
            n=1;%Number of iterations
            nEff=0;%Number of effective iterations
            absFirstIt=1;%Flag to indicate that we are in the absolute first iteration
            yF=y{s}{v}{e};
            if gpu;yF=gpuArray(yF);end
            
            for lk=1:Lk
                %SPATIAL MULTIRESOLUTION INFORMATION
                x=xGT;x(:)=0;  
                [xRes,SRes,yRes,ARes,MRes,rGridRes,kGridRes,rkGridRes]=downsampleOperators(Lk-lk,x,S,yF,A{s}{e},M,2);%Note this does not work for MB
                NXRes=size(xRes);NYRes=size(yRes);NY=size(y{s}{v}{e});
                Lteff=Lt;
                if lk~=Lk;Lteff=1;end
                downFactor=4^(Lk-lk);
                
                for lt=1:Lteff
                    contKT=contKT+1;
                    %TEMPORAL MULTIRRESOLUTION INFORMATION
                    Aeff=A{s}{e};AeffRes=ARes;Teff=TGT{s}{v};
                    if lt~=Lt && Lt>1
                        Aeff=sum(reshape(A{s}{e},[NY(1:2) 1 NS(3,s) NS(1,s)/(NS(6,s)^2) NS(4:5,s)' NS(7,s)]),4);
                        AeffRes=sum(reshape(ARes,[NYRes(1:2) 1 NS(3,s) NS(1,s)/(NS(6,s)^2) NS(4:5,s)' NS(7,s)]),4);
                        if lk==1;T=mean(reshape(T,[1 1 1 NS(3,s) NS(1,s)/(NS(6,s)^2) NS(4:5,s)' NS(7,s) 6]),4);end
                        Teff=mean(reshape(TGT{s}{v},[1 1 1 NS(3,s) NS(1,s)/(NS(6,s)^2) NS(4:5,s)' NS(7,s) 6]),4);                        
                    elseif lt>1                      
                        T=reshape(repmat(T,[1 1 1 NS(3,s) 1 1]),[1 1 1 1 prod(NS([1 3],s)/(NS(6,s)^2)) NS(4:5,s)' NS(7,s) 6]);
                    end
                    NT=size(T,5)*size(T,6)*size(T,7)*size(T,8);NC=size(S,4);
                    
                    %RECONSTRUCTION WITH KNOWN MOTION
                    x=zeros(NXRes,'like',xGT);
                    x=Xsolver(x,yRes,SRes,MRes,Teff,TFov{s},AeffRes,B{s},kGridRes,rkGridRes);
                    xGTNOTR{s}{v}{e}{contKT}=resampling(x,N)/downFactor;
                    %if contKT==2;fprintf('SNR known motion: %.2f dB\n',10*log10(sum(abs(xGT(:)).^2)/sum(abs(xGTNOTR{s}{v}{e}{contKT}(:)-xGT(:)).^2)));end
                    for f=1:2
                        if f==2;x=filtering(x,buildFilter(NXRes,'tukeyIso',1,gpu,0));end                            
                        xUp=resampling(x,N)/downFactor;
                        errXGTNOTR{s}{v}{e}(contKT,f)=gather(normm(xUp,xGT))/prod(N);%Reconstruction error
                        errFGTNOTR{s}{v}{e}(contKT,f)=gather(errorFit(xUp,yF,S,Teff,TFov{s},Aeff,B{s},kGrid,rkGrid));%Fit error with known motion                        
                    end
                    %RECONSTRUCTION FIRST ITERATION
                    [xRes,nUp]=Xsolver(xRes,yRes,SRes,MRes,T,TFov{s},AeffRes,B{s},kGridRes,rkGridRes);
                    nEff=nEff+NT*nUp/downFactor;                    
                    
                    %VARIABLES FOR XT ALGORITHM
                    EX=[];ET=[];dmin=[];
                    w=[];flagw=[];winic=[];%See directionDISORDER
                    meanT=1;
                    etRef=precomputeFactorsSincRigidTransform(kGridRes,rkGridRes,TFov{s},1,[],[],1);

                    while nEff<=nItPlan(contKT)                      
                        if debug>1;fprintf('Effective iteration: %.2f\n',nEff);end
                        
                        %ERRORS                     
                        x=resampling(xRes,N)/downFactor;
                        errX{s}{v}{e}(n)=gather(normm(x,xGT))/prod(N);
                        errF{s}{v}{e}(n)=gather(errorFit(x,yF,S,T,TFov{s},Aeff,B{s},kGrid,rkGrid));               
                        effIt{s}{v}{e}(n)=nEff;               
                        
                        %CONVERGENCE
                        %if lk==Lk && lt==Lteff;th=errFGTNOTR{s}{v}{e}(contKT,1);else th=errFGTNOTR{s}{v}{e}(contKT,2);end
                        if lk==Lk && lt==Lteff;th=errFGTNOTR{s}{v}{e}(contKT,1);else th=mean(errFGTNOTR{s}{v}{e}(contKT,:));end%Generally better criterion to illustrate
                        n=n+1;
                        if errF{s}{v}{e}(n-1)<th
                            if debug>0;fprintf('Convergence for level %d-%d reached at effective iteration %.2f\n',lk,lt,nEff);end
                            break;
                        end
                                                                       
                        %RESIDUALS
                        [et,etg]=precomputeFactorsSincRigidTransform(kGridRes,rkGridRes,T,1,1,[],1); 
                        [ry,xB]=sincRigidTransform(xRes,et,1);
                        ry=encodeDISORDER(ry,SRes,[],etRef,AeffRes,B{s},yRes);
                        nEff=nEff+NT/downFactor;
                        
                        %ENERGY COMPUTATION AND JOINT ACCEPTANCE OF MOTION UPDATE
                        EX=normm(ry,[],1:4);%if debug>1;fprintf('Error after x: %.2f\n',sum(EX));end                                                                   
                        
                        %JACOBIAN
                        if isempty(w) || any(flagw(:)==1)                            
                            Jt=sincRigidTransformGradient(xB,et,etg);
                            nEff=nEff+6*NT/downFactor;
                            xPrev=xRes;TPrev=T;ryPrev=ry;JtPrev=Jt;EXPrev=EX;
                        else
                            xRes=xPrev;T=TPrev;ry=ryPrev;Jt=JtPrev;EX=EXPrev;
                        end                                                                        

                        %GRADIENT
                        [dt,ddt]=gradientDISORDER(Jt,ry,SRes,etRef,AeffRes,B{s});

                        %MOTION UPDATE CANDIDATE
                        [Tup,w,flagw,winic]=directionDISORDER(dt,ddt,w,flagw,winic);                                               
                        NTR=size(T);
                        Tup=reshape(Tup,NTR);                  
                        Tca=T-Tup;
                        Tca=restrictTransform(Tca,[rGrid{1}(1) rGrid{1}(end);rGrid{2}(1) rGrid{2}(end);0 0]);
  
                        %RESIDUALS
                        et=precomputeFactorsSincRigidTransform(kGridRes,rkGridRes,Tca,1,[],[],1);                        
                        ry=encodeDISORDER(xRes,SRes,et,etRef,AeffRes,B{s},yRes); 
                        nEff=nEff+NT/downFactor;
                        
                            
                        %ENERGY COMPUTATION AND SEPARATED ACCEPTANCE OF MOTION UPDATE
                        ET=normm(ry,[],1:4);
                        flagw(ET(:)<=EX(:))=1;flagw(ET(:)>EX(:))=0;
                        T=resPop(dynInd(resPop(T,5:8,prod(NTR(5:8)),5),flagw==1,5,dynInd(resPop(Tca,5:8,prod(NTR(5:8)),5),flagw==1,5)),5,NTR(5:8),5:8);
                        if debug>1;fprintf('Error after T: %.2f\n',sum(min(EX(:),ET(:))));end                                               
                        
                        %SOLVE FOR X
                        ry=resPop(dynInd(resPop(ry,5:8,prod(NTR(5:8)),5),flagw==0,5,dynInd(resPop(ryPrev,5:8,prod(NTR(5:8)),5),flagw==0,5)),5,NTR(5:8),5:8);
                        [xUp,nUp]=Xsolver([],ry,SRes,MRes,T,TFov{s},AeffRes,B{s},kGridRes,rkGridRes); 
                        xRes=xRes+xUp;nEff=nEff+NT*nUp/downFactor;          
                        
                        %REMOVE AVERAGE TRANSFORM
                        if meanT
                            Tmed=multDimMea(T,5:8);
                            et=precomputeFactorsSincRigidTransform(kGridRes,rkGridRes,Tmed,1,[],[],1);
                            xRes=sincRigidTransform(xRes,et,1);
                            xRes=MRes.*xRes;
                            T=bsxfun(@minus,T,Tmed);
                        end
                    end
                end
            end     
            xEst{s}{v}{e}=gather(x);
            TEst{s}{v}{e}=gather(T);
            rEst{s}{v}{e}=errorFit(x,yF,S,T,TFov{s},Aeff,B{s},kGrid,rkGrid,1:8);%Residuals of the estimation
        end
    end
end

for s=1:size(NS,2)%Shots
    for e=1:size(recTyp,1);A{s}{e}=gather(A{s}{e});end%Encoding methods
end
xGT=gather(xGT);

