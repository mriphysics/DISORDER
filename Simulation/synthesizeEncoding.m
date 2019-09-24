function [A,B]=synthesizeEncoding(N,NS,EncMeth)

%SYNTHESIZEENCODING   Generates spectral encoding
%   [A,B]=SYNTHESIZEENCODING(NX,NS,ENCMETH)
%   generates encoding patterns included in Fig. 4 of the paper
%   * NX specifies the dimensions of the image to be reconstructed
%   * NS is a vector with first row specifying the number of shots (or 
%   of transform states), second row specifying the SENSE factors along
%   each direction and third row specifying the number of within shot
%   subdivisions
%   * ENCMETH is a cell array with each element specifying the encoding 
%   pattern (see Fig. 4 of the paper), one in 'Sequential', 'Checkered',
%   'Random-checkered', 'Random' ('Random-Stationary' and 
%   'Random-checkered-Stationary' contemplated for simulating multiple
%   dynamics but not thoroughly tested
%   ** A is a cell array of LENGTH(NS) x LENGTH(ENCMETH) containing the 
%   sampling scheme corresponding to the encoding method
%   ** B is a cell array of LENGTH(NS) x 2 containing the multiband
%   sampling scheme for each dimension
%

A=cell(size(NS,2),length(EncMeth));
B=cell(size(NS,2),2);
for S=1:size(NS,2)    
    NNS=NS(1,S);%Number of shots
    NMB=prod(abs(NS(4:5,S)));%Multiband/Slab localization factor
    NNT=NNS/NMB;%Shots per multiband/slab excitation
    NR=NS(7,S);%Number of repeats
    NSW=NS(3,S);%Within shot subdivisions
    NY=N(1:2)./(abs(NS(4:5,S))*NS(2,S))';%K-space size
    NE=NMB*prod(NY)/(NS(1,S)*NSW);%Number of echoes per shot  
    for d=1:2
        if abs(NS(3+d,S))>1
            B{S}{d}=1:N(d);
            NB=ones(1,8);
            NB(d)=N(d)/abs(NS(3+d,S));
            NB(5+d)=abs(NS(3+d,S));
            if NS(3+d,S)>0%Multiband excitation
                B{S}{d}=reshape(B{S}{d},[NB(5+d) NB(d)]);            
                B{S}{d}=B{S}{d}';            
            end
            B{S}{d}=reshape(B{S}{d},NB);%Slab excitation
        else
            B{S}{d}=[];
        end
    end    
    for E=1:length(EncMeth)        
        A{S}{E}=single(zeros([NY 1 NNT NMB NR]));
        if strcmp(EncMeth{E},'Sequential')%Sequential encoding
            for s=1:NNT
                Block=NY(2)/NNT;
                A{S}{E}(:,1+(s-1)*Block:s*Block,1,s,:,:)=1;
            end
        elseif strcmp(EncMeth{E},'Checkered')%Checkered encoding
            s=1;
            for m=1:sqrt(NNS)/abs(NS(4,S))
                BlockM=sqrt(NNS)/abs(NS(4,S));
                for n=1:sqrt(NNS)/abs(NS(5,S))
                    BlockN=sqrt(NNS)/abs(NS(5,S));
                    A{S}{E}(m:BlockM:end,n:BlockN:end,1,s,:,:)=1;
                    s=s+1;
                end
            end
        elseif strcmp(EncMeth{E},'Random')%Random encoding
            for r=1:NR
                for t=1:NMB
                    indRand=randperm(prod(NY));
                    [IR,JR]=ind2sub(NY,indRand);            
                    for s=1:NNT
                        Block=prod(NY)/NNT;
                        IRBlock=IR(1+(s-1)*Block:s*Block);
                        JRBlock=JR(1+(s-1)*Block:s*Block);
                        for k=1:length(IRBlock);A{S}{E}(IRBlock(k),JRBlock(k),1,s,t,r)=1;end
                    end
                end
            end
        elseif strcmp(EncMeth{E},'Random-Stationary')%Same randomization for different excitations
            indRand=randperm(prod(NY));
            [IR,JR]=ind2sub(NY,indRand);            
            for s=1:NNT
                Block=prod(NY)/NNT;
                IRBlock=IR(1+(s-1)*Block:s*Block);
                JRBlock=JR(1+(s-1)*Block:s*Block);
                for k=1:length(IRBlock);A{S}{E}(IRBlock(k),JRBlock(k),1,s,:,:)=1;end
            end
        elseif strcmp(EncMeth{E},'Random-checkered')%Random-checkered encoding
            NBlocks=(abs(NS(4:5,S))'.*NY)/sqrt(NNS);
            for r=1:NR
                for t=1:NMB                       
                    for o=1:NBlocks(1)
                        for p=1:NBlocks(2)
                            indRand=randperm(NNT);
                            [IR,JR]=ind2sub(sqrt(NNS)./abs(NS(4:5,S))',indRand);
                            s=1;
                            for m=1:sqrt(NNS)/abs(NS(4,S))
                                BlockM=sqrt(NNS)/abs(NS(4,S));
                                for n=1:sqrt(NNS)/abs(NS(5,S))
                                    BlockN=sqrt(NNS)/abs(NS(5,S));
                                    A{S}{E}(IR(s)+(o-1)*BlockM,JR(s)+(p-1)*BlockN,1,s,t,r)=1;
                                    s=s+1;
                                end
                            end
                        end
                    end
                end
            end
        elseif strcmp(EncMeth{E},'Random-checkered-Stationary')%Same randomization for different excitations
            NBlocks=(abs(NS(4:5,S))'.*NY)/sqrt(NNS);                     
            for o=1:NBlocks(1)
                for p=1:NBlocks(2)
                    indRand=randperm(NNT);
                    [IR,JR]=ind2sub(sqrt(NNS)./abs(NS(4:5,S))',indRand);
                    s=1;
                    for m=1:sqrt(NNS)/abs(NS(4,S))
                        BlockM=sqrt(NNS)/abs(NS(4,S));
                        for n=1:sqrt(NNS)/abs(NS(5,S))
                            BlockN=sqrt(NNS)/abs(NS(5,S));
                            A{S}{E}(IR(s)+(o-1)*BlockM,JR(s)+(p-1)*BlockN,1,s,:,:)=1;
                            s=s+1;
                        end
                    end
                end
            end
        else
            error('Undefined Encoding method');
        end
        for m=1:2;A{S}{E}=ifftshift(A{S}{E},m);end%DC at first element
        A{S}{E}=reshape(A{S}{E},[NY 1 1 NNS NR]);
        NNSR=round(NNS/(NS(6,S))^2);
        NNTR=round(NNT/(NS(6,S))^2);
        A{S}{E}=A{S}{E}(:,:,:,:,1:NNSR,:);

        C=repmat(A{S}{E},[1 1 1 NSW 1 1]);
        C(:)=0;      
        for r=1:NR
            for s=1:NNSR
                cont=0;
                contE=1;
                for m=1:NY(1)
                    for n=1:NY(2)                    
                        if A{S}{E}(m,n,1,1,s,r)==1
                            C(m,n,1,contE,s,r)=1;
                            cont=cont+1;
                            if cont==NE
                                contE=contE+1;
                                cont=0;
                            end
                        end
                    end
                end
            end
        end        
        A{S}{E}=reshape(C,[NY 1 1 NNTR*NSW abs(NS(4:5,S))' NR]);
    end
end
