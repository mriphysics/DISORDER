function [y,xGTNO,xNOTR,errFGTNO,errXGTNO,no]=synthesizeY(xGT,TGT,TFov,S,A,B,M,kGrid,rkGrid,NS,sigma,no)

%SYNTHESIZEY   Synthesizes motion corrupted data
%   [Y,XGTNO,XGTNOTR,ERRFGTNO,ERRXGTNO,NO]=SYNTHESIZEY(XGT,TGT,TFOV,S,A,B,M,KGRID,RKGRID,NS,{SIGMA},{NO})
%   synthesizes motion corrupted data according to the forward model in 
%   Fig. 1 of the TCI paper,
%   * XGT is the ground truth image
%   * TGT is a cell array of cell array of #TestedShots x #MotionLevels 
%   ground truth transforms
%   * TFOV is a cell array of cell array of #TestedShots x #MotionLevels 
%   FOV transforms (for multi-oriented data)
%   * S is the coil array sensitivity map
%   * A is a cell array of #TestedShots x #EncodingMethods containing the 
%   sampling scheme
%   * B is a cell array of #TestedShots x 2 containing the MB sampling 
%   scheme for each dimension
%   * M is a mask
%   * KGRID is a 1x3 cell array with the spectral coordinates along each 
%   axis
%   * RKGRID is a 2x3 cell array where the first dimension of the cell
%   indexes the shearing orientation and the second indexes the 
%   spatio-spectral axes pair used for that particular shearing (given by 
%   [1 3 2;2 1 3]). See ec. 4 of the paper for further details.
%   * NS is a vector with first row specifying the number of shots (or 
%   of transform states), second row specifying the SENSE factors along
%   each direction and third row specifying the number of within shot
%   subdivisionsplugNoise(repmat(S,[1 1 1 1 NSS NMB' max(NR)]))
%   * {SIGMA} is the noise level to be added, it defaults to 0
%   * {NO} is a previously synthesized noise sample
%   ** Y is the synthesized data
%   ** XGTNO is the ground-truth reconstruction in the presence of noise
%   ** XNOTR is the reconstruction in the presence of noise without motion
%   estimation
%   ** ERRFGTNO is the fitting error of the ground-truth reconstruction in
%   the presence of noise
%   ** ERRXGTNO is the error of the ground-truth reconstruction in the 
%   presence of noise
%   ** NO is the estimated noise
%

NMB=max(abs(NS(4:5,:)),[],2);
NSS=max(NS(1,:));

if nargin<11 || isempty(sigma);sigma=0;end
if nargin<12 || isempty(no) || size(no,4)~=size(S,4);no=sigma*plugNoise(repmat(S,[1 1 1 1 NSS NMB' max(NS(7,:))]));end

N=size(xGT);N(end+1:3)=1;%Image size
gpu=isa(xGT,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end

for s=1:length(TGT)
    NY=size(A{s}{1});NY(end+1:8)=1;
    NNR=NS(7,s);
    %Synthesis in the absence of motion and acceleration  
    yMF=bsxfun(@plus,encodeDISORDER(xGT,S,[],[],A{s}{1},B{s}),dynInd(no,{1:NY(1),1:NY(2),1:NY(5),1:abs(NS(4,s)),1:abs(NS(5,s)),1:NS(7,s)},[1:2 5:8]));
    %Reconstruction in the absence of motion
    x=zeros(N,'like',xGT);
    xGTNO{s}{1}=Xsolver(x,yMF,S,M,[],[],A{s}{1},B{s});
    
    yAlt=bsxfun(@plus,encodeDISORDER(xGT,S,[],[],A{s}{1},B{s}),dynInd(no,{1:NY(1),1:NY(2),1:NY(5),1:abs(NS(4,s)),1:abs(NS(5,s)),1:NS(7,s)},[1:2 5:8])*sqrt(NY(8)));
    %Reconstruction in the absence of motion-single NSA
    x=zeros(N,'like',xGT);
    xGTNO{s}{2}=Xsolver(x,yAlt,S,M,[],[],A{s}{1},B{s});
    
    yAlt=bsxfun(@plus,encodeDISORDER(xGT,S,[],[],A{s}{1},B{s}),dynInd(no,{1:NY(1),1:NY(2),1:NY(5),1:abs(NS(4,s)),1:abs(NS(5,s)),1:NS(7,s)},[1:2 5:8])*sqrt(NY(8)*NS(1,s)));
    %Reconstruction in the absence of motion-single NSA
    x=zeros(N,'like',xGT);
    xGTNO{s}{3}=Xsolver(x,yAlt,S,M,[],[],A{s}{1},B{s});
    
    
    %Fit error in the absence of motion
    errFGTNO{s}=gather(errorFit(xGTNO{s}{1},yMF,S,[],[],A{s}{1},B{s}));
    errXGTNO{s}=gather(normm(xGTNO{s}{1},xGT));
    for l=1:3;xGTNO{s}{l}=gather(xGTNO{s}{l});end
    
    etRef=precomputeFactorsSincRigidTransform(kGrid,rkGrid,TFov{s},1,0,[],1);
    for v=1:length(TGT{s})        
        et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,TGT{s}{v},1,0,[],1);
        for e=1:length(A{s})
            %SYNTHESIS
            NY=size(A{s}{e});NY(end+1:8)=1;      
            x=encodeDISORDER(xGT,S,et,etRef,A{s}{e},B{s});
            x=sum(x,5);
            y{s}{v}{e}=gather(x+no(1:NY(1),1:NY(2),:,:,1:NY(5),1:NY(6),1:NY(7),1:NY(8)));
                        
            %RECONSTRUCTION WITH ZERO MOTION
            xx=zeros(N,'like',xGT);
            xNOTR{s}{v}{e}=gather(Xsolver(xx,y{s}{v}{e},S,M,[],TFov{s},A{s}{e},B{s},kGrid,rkGrid));                          
        end        
    end
end
no=gather(no);

