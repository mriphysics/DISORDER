function [k,NP,K,kGrid,ShotNumber,EchoNumber]=samplingDISORDER(k,U,SegmentOrder,SampleOrder,K)

%SAMPLINGDISORDER is a template to compute the motion-sensitive sampling 
%orders in [1] L Cordero-Grande et al, "Motion corrected magnetic resonance 
%imaging %with DISORDER: Distributed and Incoherent Sample Orders for 
%Reconstruction Deblurring using Encoding Redundancy", 2019.
%   K=SAMPLINGDISORDER(K,{U},{SEGMENTORDER},{SAMPLEORDER},K)
%   * K is the set of phase encoding (PE) sampling points in any given 
%   order, with the two phase encodes along the columns and the samples 
%   along the rows
%   * {U} is the tiling size. It defaults to round(sqrt(NP)) X 
%   round(sqrt(NP))
%   * {SEGMENTORDER} is the ordering used to acquire the segments. Possible 
%   orderings are 'Sequential','Checkered','RandomCheckered' and 'Random'. 
%   It defaults to 'RandomCheckered'
%   * {SAMPLEORDER} is the ordering used to acquire the samples within each
%   segment. Possible orderings are 'ZigZag' and 'AlternatingZigZag'. It 
%   defaults to AlternatingZigZag
%   * {K} serves to set the spectral limits, mainly for visualization
%   purposes, otherwise it is obtained from the sampling points KIN
%   ** K is the reordered trajectory
%   ** NP is the number of profiles
%   ** K are the spectral limits
%   ** KGRID is the spectral grid as a cell of size 1x2 with each element
%   including the grid along a given dimension (respectively as a row and a column 
%   vector)
%   ** SHOTNUMBER is the shot index for each sample
%   ** ECHONUMBER is the echo index for each sample
%   

NP=size(k,1);%Number of profiles
if nargin<2 || isempty(U);U=round(sqrt(NP))*ones(1,2);end
if nargin<3 || isempty(SegmentOrder);SegmentOrder='RandomCheckered';end
if nargin<4 || isempty(SampleOrder);SampleOrder='AlternatingZigZag';end
if nargin<5 || isempty(K);K=max(k,[],1)-min(k,[],1)+1;end
fprintf('Spectral dimensions: %d %d\n',K(1),K(2));

S=prod(U);%Number of shots
fprintf('Number of shots: %d\n',S);

Uor=U;SegmentOrderOr=SegmentOrder;
if strcmp(SegmentOrder,'Sequential');SegmentOrder='Checkered';U=[1 1];k=flip(k,2);K=flip(K);end
S=prod(U);%Effective number of shots

assert(strcmp(SampleOrder,'AlternatingZigZag') || mod(NP,S)==0,'Number of profiles (%d) not a multiple of number of shots (%d)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRECOMPUTATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kGrid=generateGrid(K,0,K,ceil((K+1)/2),0);

if strcmp(SegmentOrder,'Random') || strcmp(SegmentOrder,'RandomCheckered') || strcmp(SegmentOrder,'Checkered')
    %%%%%%%%%%%%%%%%%%%%%%%GENERATING THE TILE STRUCTURE%%%%%%%%%%%%%%%%%%%%%%%%
    Ul=cell(1,2);
    UperK=zeros(1,2);
    centU=zeros(1,2);
    for d=1:2
        Ul{d}=kGrid{d}(mod(kGrid{d}-ceil(U(d)/2),U(d))==0);
        Ul{d}=Ul{d}(:);
        Ul{d}=cat(1,Ul{d}(1)-diff(Ul{d}(1:2)),Ul{d});
        Ul{d}(:,2)=Ul{d}+U(d)-1;
        Ul{d}=Ul{d}';
        UperK(d)=size(Ul{d},2);
        centU(d)=find(Ul{d}(2,:)>0,1,'first');
    end
 
    %%%%%%%%%%%%%%%%%%%%%%DETECTING TILES INSIDE THE SHUTTER%%%%%%%%%%%%%%%%%%%
    mk=min(k,[],1)-1;
    ElSh=zeros(prod(K),1);
    ElSh(sub2indV(K,bsxfun(@minus,k,mk)))=1;%Eliptical shutter
    NPBlock=zeros(UperK);
    lim=zeros(2,2);
    for m=1:UperK(1)
         for n=1:UperK(2)
             po=[m n];         
             for d=1:2;lim(:,d)=Ul{d}(:,po(d));end
             if lim(1,1)>=kGrid{1}(1) && lim(2,1)<=kGrid{1}(end) && lim(1,2)>=kGrid{2}(1) && lim(2,2)<=kGrid{2}(end)
                 corn=[lim(1,1) lim(1,2); lim(2,1) lim(1,2); lim(1,1) lim(2,2); lim(2,1) lim(2,2)];         
                 NPBlock(m,n)=sum(ElSh(sub2indV(K,bsxfun(@minus,corn,mk))));
             end
         end
    end
    NPBlock=NPBlock==4;%All corners inside

    %%%%%%%%%%%%%%%LABEL POINTS INSIDE CHECKERBOARD BLOCKS%%%%%%%%%%%
    posAcc=ind2subV(UperK,1:prod(UperK));
    posAcc=posAcc(NPBlock(:),:);    
    assignProf=zeros(1,NP);
    for m=1:size(posAcc,1)
        for d=1:2;lim(:,d)=Ul{d}(:,posAcc(m,d));end
        assignProf(all(bsxfun(@ge,k,lim(1,:)) & bsxfun(@le,k,lim(2,:)),2))=1;
    end    
    NPOu=sum(~assignProf);
    kk=zeros(NP,2);

    %%%%%%%%%%%%POINTS CORRESPONDING TO BLOCKS INSIDE THE SHUTTER%%%%%%%%%%%%%%
    cont=1;
    for pp=1:size(posAcc,1)
        kk(cont:cont+S-1,:)=bsxfun(@plus,ind2subV(U,1:S),[Ul{1}(1,posAcc(pp,1)) Ul{2}(1,posAcc(pp,2))]-1);
        cont=cont+S;
    end
    %%%%BLOCKS AT THE BOUNDARY OF THE SHUTTER (ANGULAR TRAVERSAL FOR SHOTS)%%%%
    KKOut=single(k(~assignProf,:));
    KKOutAng=atan2(KKOut(:,2),KKOut(:,1));%Angular k-coordinates
    KKOutAng=wrapToPi(-KKOutAng);
    KKOutRad=sqrt(sum(KKOut.^2,2));%Radial k-coordinates
    corF=1e-5;%This gives good numerical results. It just bends a tiny bit the radial coordinates
    [~,indRsort]=sort(KKOutAng(:)+corF*KKOutRad(:));%Angular sorting of k-coordinates
    kk(cont:cont+NPOu-1,:)=KKOut(indRsort,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%CHECKER REORDERING%%%%%%%%%%%%%%%%%%%%%%%%
    Srand=electrostaticRepulsionDISORDER(U);
    %Srand=1:S;%Note that uncommenting this line would disable the within tile distribution of samples
    for e=1:S:S*floor(NP/S)
        vE=e:e+S-1;
        if strcmp(SegmentOrder,'RandomCheckered');Srand=randperm(S);end
        kk(vE(Srand),:)=kk(vE,:);
    end
    if strcmp(SegmentOrder,'Random');kk=kk(randperm(NP),:);end        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%ECHO ORDERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for s=1:S
        vS=s:S:NP;
        kke=kk(vS,:);
        indBlKKE=zeros(size(kke,1),2);
        for l=1:size(kke,1)
            ent=0;
            for m=1:UperK(1)
                if kke(l,1)>=Ul{1}(1,m) && kke(l,1)<=Ul{1}(2,m)
                    for n=1:UperK(2)
                        if kke(l,2)>=Ul{2}(1,n) && kke(l,2)<=Ul{2}(2,n)
                            indBlKKE(l,:)=[m n];
                            if mod(m,2)==0
                                kke(l,2)=-kke(l,2);
                                indBlKKE(l,2)=UperK(2)+1-n;
                                ent=1;
                                break;
                            end
                        end
                    end
                    if ent;break;end
                end
            end
        end
        if any(indBlKKE(:))==0;error('Sample not assigned to a block');end
        or1=1e9*indBlKKE(:,1);
        or2=1e6*indBlKKE(:,2);
        if strcmp(SegmentOrderOr,'Sequential') && strcmp(SampleOrder,'ZigZag');or2=or2(randperm(length(or2)));end
        or3=1e3*kke(:,1);
        or4=kke(:,2);
        [~,indEchoSort]=sort(or1+or2+or3+or4);
        if mod(s,2)==0 && strcmp(SampleOrder,'AlternatingZigZag');indEchoSort=flip(indEchoSort);end%For Eddy currents in non-shot sequences
        if ~strcmp(SegmentOrder,'Random');kk(vS,:)=kk(vS(indEchoSort),:);end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%FINAL TRAJECTORY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=zeros(size(kk));
    cont=1;
    for s=1:S
        vS=s:S:NP;
        k(cont:cont+length(vS)-1,:)=kk(vS,:);
        cont=cont+length(vS);
    end
end

if strcmp(SegmentOrderOr,'Sequential')
    k=flip(k,2);K=flip(K);
    kGrid=generateGrid(K,0,K,ceil((K+1)/2),0);
    U=Uor;
    S=prod(U);
end

if strcmp(SegmentOrderOr,'Sequential') && strcmp(SampleOrder,'ZigZag')
    kk=k;
    NE=NP/S;    
    for s=1:S
        vS=(s-1)*NE+1:s*NE;
        kke=kk(vS,:);
        or1=1e3*kke(:,1);
        or2=kke(:,2);   
        or2(mod(kke(:,1),2)==0)=-or2(mod(kke(:,1),2)==0);
        [~,indEchoSort]=sort(or1+or2);
        kk(vS,:)=kk(vS(indEchoSort),:);
    end 
    k=kk;
end

if nargout>1
    S=prod(Uor);
    EchoNumber=zeros(K);
    ShotNumber=zeros(K);
    cont=1;
    for s=1:S
        vS=s:S:NP;
        for e=1:length(vS)
            EchoNumber(k(cont,1)-kGrid{1}(1)+1,k(cont,2)-kGrid{2}(1)+1)=e;
            ShotNumber(k(cont,1)-kGrid{1}(1)+1,k(cont,2)-kGrid{2}(1)+1)=s;
            cont=cont+1;
        end
    end
end


