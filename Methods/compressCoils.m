function [S,y,C,D,U]=compressCoils(S,perc,y,useFold,useInvMot)

% COMPRESSCOILS performs channel compression based on [1] M Buehrer, 
%   KP Pruessmann, P Boesiger, S Kozerke, "Array Compression for MRI With 
%   Large Coil Arrays," Magn Reson Med, 57:1131-1139, 2007
%   [S,Y,C,D,U]=COMPRESSCOILS(S,PERC,{Y},{USEFOLD},{USEINVMOT}) 
%   * S is the original channel information
%   * PERC is the energy rate to be kept if strictly lower than one /
%   number of components to be kept if bigger or equal than one. If empty
%   no compression is performed
%   * {Y} is the original data sampled across channels
%   * {USEFOLD} is a flag that denotes whether to use the folding 
%   structure when compressing. Defaults to 1 if {Y} is present and to {0} 
%   otherwise
%   * {USEINVMOT} is a flag that indicates to use the inverse coils for 
%   motion estimation, deprecated
%   ** S is the compressed channel information
%   ** Y is the compressed data
%   ** C are the number of components to use for reconstruction (first 
%   element) and motion estimation (last component)
%   ** D are the singular values of the decomposition
%   ** U are the singular vectors of the decomposition
%

if nargin<3;y=[];end
if nargin<4 || isempty(useFold);useFold=1;end
if nargin<5 || isempty(useInvMot);useInvMot=0;end

gpu=isa(S,'gpuArray');gpuIn=single(gpuDeviceCount && ~blockGPU);

reg=1e-9;

%Initialize
if ~isempty(y);NY=size(y);else NY=size(S);end;NY(end+1:14)=1;
NS=size(S);NS(end+1:4)=1;

assert(ndims(S)<=4,'Dimensionality of coil profiles (%d) should not be bigger than 4',ndims(S));

C=NS(4)*ones(1,length(perc));
if ~isempty(perc)
    %Compute P
    if useFold==1 && any(NS(1:3)~=NY(1:3))
        x=cell(1,3);vM=cell(1,3);     
        for m=1:3      
            x{m}=single(ones([NS(m) 1]));
            x{m}=fold(x{m},1,NS(m),NY(m));
            x{m}=ifold(x{m},1,NS(m),NY(m));
            vM{m}=unique(x{m});
        end
        P=single(zeros(NY(4),NY(4)));
        if gpuIn;P=gpuArray(P);end
        for m=1:length(vM{1})
            for n=1:length(vM{2})
                for o=1:length(vM{3})
                    SF=dynInd(S,{x{1}==vM{1}(m),x{2}==vM{2}(n),x{3}==vM{3}(o)},1:3); 
                    if gpuIn;SF=gpuArray(SF);end
                    NSF=size(SF);NSF(end+1:4)=1;
                    SF=resPop(SF,[1 2 3 4],[NSF(1)/vM{1}(m) vM{1}(m) NSF(2)/vM{2}(n) vM{2}(n) NSF(3)/vM{3}(o) vM{3}(o) NSF(4)],[5 2 6 3 7 4 1]);
                    NSF=size(SF);NSF(end+1:7)=1;                
                    SF=reshape(SF,[NSF(1) prod(NSF(2:4)) prod(NSF(5:7))]);
                    BlSz=1e5;
                    for p=1:BlSz:size(SF,3);vn=p:min(p+BlSz-1,size(SF,3));
                        SFaux=dynInd(SF,vn,3);                                       
                        SFH=matfun(@ctranspose,SFaux);
                        Sbis=matfun(@mtimes,SFH,SFaux);
                        Sbis=matfun(@mldivide,Sbis,SFH);
                        Sbis=matfun(@mtimes,SFaux,Sbis);
                        P=P+sum(Sbis,3);
                    end
                 end                 
             end
        end
        Sbis=[];SF=[];SFH=[];
    else
        SH=conj(S);
        Normal=(sum(SH.*S,4)+reg).^(-1);
        Saux=bsxfun(@times,S,Normal);Normal=[];
        perm=1:14;perm([4 5])=[5 4];
        SH=permute(SH,perm);
        P=single(zeros(NY(4),NY(4)));
        for l=1:NS(3)
            Sbis=bsxfun(@times,dynInd(Saux,l,3),dynInd(SH,l,3));    
            P=P+shiftdim(multDimSum(Sbis,1:2),3);
        end
        Saux=[];SH=[];Sbis=[];
    end
    P=gather(P);
    %Compute F
    [U,D]=svd(P);   
    D=diag(abs(D));
    %Thresholding 
    if useInvMot
        if perc(1)<1;C(1)=find(cumsum(D)/sum(D)>=perc(1),1);else C(1)=perc(1);end%Thresholding/Preservation
        vM=1:C(1);
        if length(perc)>=2
            if perc(2)<1;C(2)=find(flip(cumsum(flip(D)))/sum(D)<1-perc(2),1);else C(2)=perc(2);end
            vM=horzcat(vM,max(C(1)+1,C(2)):NS(4));
            C(2)=NS(4)-C(2)+1;
            if length(perc)==3
                if perc(1)<1;C(3)=find(cumsum(D)/sum(D)>=perc(3),1);else C(3)=perc(3);end%Thresholding/Preservation
            end        
        end
    else
        for n=1:length(perc)
            if perc(n)<1;C(n)=find(cumsum(D)/sum(D)>=perc(n),1);else C(n)=perc(n);end%Thresholding/Preservation
        end
        vM=1:max(C);
    end
    %fprintf('Number of components for final reconstruction / motion estimation / reconstruction for motion estimation:%s\n',sprintf(' %d',C));    
    
    A=dynInd(ctranspose(U),vM,1);%Compressing matrix
    if gpu;A=gpuArray(A);end
    S=aplGPU(A,S,4);%Compress coils   
    if ~isempty(y);y=aplGPU(A,y,4);end%Compress data
    if gpu;U=gpuArray(U);D=gpuArray(D);end
end
