function fig0708(id,pr)

%FIG0708  Generates Figs. 7, and 8 of the manuscript "Motion corrected
%magnetic resonance imaging with DISORDER: Distributed and Incoherent 
%Sample Orders for Reconstruction Demixing using Encoding Redundancy," L 
%Cordero-Grande, G Ferrazzi, RPAG Teixeira, AN Price, J O'Muircheartaigh, 
%and JV Hajnal, arXiv:1910.00540, 2019
%   FIG0708({ID},{PR})
%   * {ID} indicates the figure to generate, either '07' or '08', if
%   empty (default), it generates all figures
%   * {PR} calls the reconstruction and visualization (0), only the
%   reconstruction (1) or only the visualization (2). It defaults to 1
%

pathInpData='../DISORDERData';%Path with input data, modify if required
pathOutData7='../DISORDERData/Results7';%Path with output data, modify if required
pathOutData8A='../DISORDERData/Results8-A';%Path with output data, modify if required
pathOutData8B='../DISORDERData/Results8-B';%Path with output data, modify if required
if ~exist(pathOutData7,'dir');mkdir(pathOutData7);end
if ~exist(pathOutData8A,'dir');mkdir(pathOutData8A);end
if ~exist(pathOutData8B,'dir');mkdir(pathOutData8B);end

if nargin<1 || isempty(id);id={'07','08'};elseif ~iscell(id);id={id};end
if nargin<2 || isempty(pr);pr=1;end

addpath(genpath('..'));

namF{1}={'GT','Q1','Q2','Q3'};
namF{2}={'MPRAGE','TSE','FLAIR','SPGR','BSSFP'};

for ii=1:length(id)
    clearvars -except id pr pathInpData pathOutData7 pathOutData8A pathOutData8B ii namF    
    gpu=(gpuDeviceCount>0 && ~blockGPU);%0->Use CPU / 1->Use GPU
    
    %EXPERIMENT ID
    figN=id{ii};%From '07' to '08'
    if strcmp(figN,'07');nam=namF{1};elseif strcmp(figN,'08') nam=namF{2};else fprintf('Figure %s not contemplated\n',id{ii});return;end
    if pr~=2    
        for n=1:length(nam)
            tsta=tic;
            load(sprintf('%s/%s.mat',pathInpData,nam{n}),'recInp');
            tend=toc(tsta);fprintf('\n\nTime loading case %s-%s: %.3f s\n',figN,nam{n},tend);
            recInp.Alg.WriteSnapshots=0;%Set to 1 to enable visualization
            
            %FOR CONSISTENCY TESTS, SHOULD GENERATE ALL 0 NIFTIS
            %recInp.Alg.parXT.exploreMemory=1;
            
            %FOR QUICK INSPECTION-ALMOST NO CORRECTION BUT SHOULD GENERATE
            %SOME NIFTI VOLUMES
            %recInp.Alg.parXT.traLimX=recInp.Alg.parXT.traLimX*100;
            %recInp.Alg.parXT.traLimXT=recInp.Alg.parXT.traLimXT*100;
            %recInp.Alg.parXT.NWend=min(recInp.Alg.parXT.NWend,2);
            if strcmp(figN,'07')
                recInp.Names.pathOu=pathOutData7;                
                recInp.Alg.AlignedRec=4;%Within shot
                recInp.Alg.parXT.computeCSRecon=0;%Not computing regularized reconstruction
                solveXT(recInp);
            else
                recInp.Names.pathOu=pathOutData8A;
                recInp.Alg.AlignedRec=1;%Non-robust
                recInp.Alg.parXT.computeCSRecon=0;%Not computing regularized reconstruction
                solveXT(recInp);
                recInp.Names.pathOu=pathOutData8B;
                recInp.Alg.AlignedRec=2;%Robust
                recInp.Alg.parXT.computeCSRecon=1;%Computing regularized reconstruction
                solveXT(recInp);
            end
            %return
        end
    end
    if pr~=1
        if strcmp(figN,'07')
            suff={'Aq','Di'};            
            for s=1:length(nam)
                fileName=sprintf('%s/An-Ve/%s',pathOutData7,nam{s});
                x=readNII(fileName,suff,gpu);
                y{1}=abs(x{1});
                y{2}=abs(dynInd(x{2},size(x{2},4),4));
                y{3}=abs(dynInd(x{2},1,4));
                x=[];
                z{s}=cat(5,y{:});y=[];%Dimension 5 is different levels of correction
            end
            w=cat(4,z{:});z=[];%Dimension 4 is different corruptions
            NW=size(w);

            %REGISTRATION
            w=resSub(w,4:5);
            W=w;W=W(:,:,:,1);W(:)=1;W([1:floor(NW(1)/12) 7*ceil(NW(1)/12):NW(1)],:,:)=0;

            pyr=[32 16 8 4 2 1];
            w=gather(abs(groupwiseVolumeRegistration(w,W,[],[],pyr)));
            w=reshape(w,NW);

            perm=[3 2 1 4 5];
            w=permute(w,perm);
            w=abs(flip(permute(dynInd(w,1:floor((1-1/3)*NW(1)),3),[2 1 3 4 5]),2));
            w=dynInd(w,round(size(w,3)/2),3);
            w=dynInd(w,{4:117,6:142},1:2);

            e=abs(dynInd(w,2:4,4)-dynInd(w,1,4));

            limInt=[0.999 0.995];
            z=sort(w(:));
            z=z(round(length(z)*limInt(1)));
            w=min(w,z);
            w=double(gather(w/z));
            w=flip(flip(w,1),2);

            z=sort(e(:));
            z=z(round(length(z)*limInt(2)));
            e=min(e,z);
            e=double(gather(e/z));
            e=flip(flip(e,1),2);

            w=cat(4,dynInd(w,1:2,4),dynInd(e,1,4),dynInd(w,3,4),dynInd(e,2,4),dynInd(w,4,4),dynInd(e,3,4));
            w=permute(w,[2 1 3 4 5]);
            w=flip(w,1);
            NW=size(w);
            
            w=permute(w,[1 5 2 4 3]);
            w=reshape(w,NW(1)*NW(5),[]);
            pause(1);figure;pause(1)            
            imshow(w,[]);%,'Border','tight')
            set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])
            title('See caption of Fig. 7 in the manuscript','FontSize',20)
        else            
            suff={'Aq','Di','Re'};
            limInt=[1 0.98 0.995 0.9925 0.95];
            for s=1:length(nam)
                fileName=sprintf('%s/An-Ve/%s',pathOutData8A,nam{s});
                x=readNII(fileName,suff(1:2),gpu);
                fileName=sprintf('%s/An-Ve/%s',pathOutData8B,nam{s});
                xaux=readNII(fileName,suff(3),gpu);
                x=[x xaux];xaux=[];                
                for l=1:3
                   y{l}=dynInd(x{l},1,4);
                   perm=[3 2 1];
                   y{l}=permute(y{l},perm);
                   y{l}=abs(flip(permute(dynInd(y{l},1:floor((1-1/3)*size(y{l},3)),3),[2 1 3]),2));
                   y{l}=dynInd(y{l},round(size(y{l},3)/2),3)';      
                   y{l}=y{l}(23:222,21:170);
                   z=sort(y{l}(:));
                   z=z(round(length(z)*limInt(s)));
                   y{l}=min(y{l},z);
                   y{l}=double(gather(y{l}/z));                 
                end
                w{s}=cat(5,y{:});
            end
            z=cat(4,w{:});
            NZ=size(z);
            z=permute(z,[1 5 2 4 3]);
            z=reshape(z,NZ(1)*NZ(5),[]);
            pause(1);figure;pause(1)
            imshow(z,[])%,'Border','tight')
            set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])
            title('See caption of Fig. 8 in the manuscript','FontSize',20)
        end            
    end
end
