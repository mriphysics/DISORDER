function fig0506(id,quick,pr)

%FIG050607  Generates Figs. 5 and 6 of the manuscript "Motion corrected
%magnetic resonance imaging with DISORDER: Distributed and Incoherent 
%Sample Orders for Reconstruction Demixing using Encoding Redundancy," L 
%Cordero-Grande, G Ferrazzi, RPAG Teixeira, AN Price, J O'Muircheartaigh, 
%and JV Hajnal, arXiv:1910.00540, 2019
%   FIG0506({ID},{ACCEL},{PR})
%   * {ID} indicates the type of study, either 'Encoding', 
%   'Multiresolution' or 'Acceleration', if empty (default), it generates 
%   all studies
%   * {QUICK} if activated, performs a quick run for simplified versions of
%   the figures in the manuscript involving 5x less iterations, not 
%   including the right-most column and bottom-most row of the convergence
%   plots. If 0 it runs the tests in the manuscript
%   * {PR} flags the generation of a new set of motion parameters / noise
%   realizations and performs the simulations (0, default), uses previously 
%   generated parameters / realizations to perform the simulations (1) or 
%   plots results of previous simulation (2)
%

pathInpData='../DISORDERData';%Path with input data, modify if required
pathOutData='../DISORDERData/Results5-6';%Path with output data, modify if required
if ~exist(pathOutData,'dir');mkdir(pathOutData);end

if nargin<1 || isempty(id);id={'Encoding','Multiresolution','Acceleration'};elseif ~iscell(id);id={id};end
if nargin<2 || isempty(quick);quick=1;end%SET THIS PARAMETER TO 0 TO FULLY REPRODUCE THE EXPERIMENTS IN THE MANUSCRIPT
if nargin<3 || isempty(pr);pr=0;end

for n=1:length(id)
    if strcmp(id{n},'Encoding');id{n}='05';
    elseif strcmp(id{n},'Multiresolution');id{n}='06';
    elseif strcmp(id{n},'Acceleration');id{n}='07';
    end
end

addpath(genpath('..'));

for ii=1:length(id)
    clearvars -except id quick pr pathInpData pathOutData ii

    %EXPERIMENT ID
    expN=id{ii};%From '05' to '07'
    if ~(strcmp(expN,'05') || strcmp(expN,'06') || strcmp(expN,'07'));fprintf('Figure %s not contemplated\n',id{ii});return;end

    %GPU/DEBUG INFO
    gpu=(gpuDeviceCount>0 && ~blockGPU);%0->Use CPU / 1->Use GPU
    if gpu;dev=gpuDevice;end
    if gpu;gpuF=2;else gpuF=0;end
    debug=2;%0/1/2 increasing the amount of debug info provided 

    if ismember(pr,0:1)
        %LOAD SYNTHETIC DATA
        % - Ground truth image xGT of size 128x128
        % - Sensitivity maps S of size 128x128x1x32
        load(sprintf('%s/xGT',pathInpData),'xGT','S');
        if gpu;xGT=gpuArray(xGT);S=gpuArray(S);end

        %EXPERIMENT PARAMETERS
        [theta,NS,recTyp]=generateParametersExp(expN,quick);

        %GRID GENERATION
        N=size(xGT);N(end+1:3)=1;%Image size
        cent=ceil((N+1)/2);
        [rGrid,kGrid,rkGrid,kkGrid]=generateTransformGrids(N,gpu,N,cent,1);

        %LOAD/SYNTHESIZE MOTION/ENCODING/NOISE
        if pr
            load(sprintf('%s/Fig%s',pathOutData,expN),'TGT','TFov','no','A','B');
            if gpu
                for s=1:size(NS,2)%Shots
                    for b=1:2;B{s}{b}=gpuArray(B{s}{b});end               
                    for e=1:size(recTyp,1);A{s}{e}=gpuArray(A{s}{e});end%Encoding methods
                end
                no=gpuArray(no);
            end
        else
            [A,B]=synthesizeEncoding(N,NS,recTyp(:,1));
            [TGT,TFov]=synthesizeT(NS,theta);
            no=[];
        end

        %ISOTROPIC RESOLUTION
        xGT=filtering(xGT,buildFilter(N,'tukeyIso',1,gpu,0));

        %COIL ARRAY COMPRESSION TO ACCELERATE EXPERIMENTS
        perc=0.99;
        S=compressCoils(S,perc);

        %SPATIAL MASK (disabled)
        M=ones(N,'like',real(xGT));
        xGT=xGT.*M;

        %NORMALIZE
        xGT=xGT/norm(xGT(:));

        %DATA SYNTHESIS
        %snrdB=[30 20 10];%High / Acceptable / Poor SNR
        snrdB=30;
        sigma=SNRToLevels(snrdB,xGT,S);
        [y,xGTNO,xNOTR,errFGTNO,errXGTNO,no]=synthesizeY(xGT,TGT,TFov,S,A,B,M,kGrid,rkGrid,NS,sigma,no);
    
        %ITERATIONS
        if strcmp(expN,'05');nExtern=20000;
        elseif strcmp(expN,'06');nExtern=200000;
        elseif strcmp(expN,'07');nExtern=200000;
        end%Maximum number of iterations of the joint method
        if quick;nExtern=ceil(nExtern/5);end
        
        %SOLVER CALL
        pyramidDISORDER;
        
        %WRITE DATA
        save(sprintf('%s/Exp%s',pathOutData,expN),'xNOTR','xEst','xGT','xGTNO','xGTNOTR','TEst','TGT','TFov','rEst','A','B','no','errFGTNO','errXGTNO','errFGTNOTR','errXGTNOTR','effIt','errX','errF','theta','NS','recTyp','nExtern','snrdB');
    else
        load(sprintf('%s/Exp%s',pathOutData,expN),'xNOTR','xEst','xGT','xGTNO','xGTNOTR','TEst','TGT','TFov','rEst','A','B','no','errFGTNO','errXGTNO','errFGTNOTR','errXGTNOTR','effIt','errX','errF','theta','NS','recTyp','nExtern','snrdB');
    end
        
    if strcmp(expN,'07');dfv=1;else dfv=0;end

    E=size(recTyp,1);
    V=size(theta,2)-dfv;
    S=size(NS,2);

    cont=1;

    FontSizeA=14;%10;
    FontSizeB=22;%18;
    FontSizeC=32;%28;
    FontSizeD=10;
    FontSizeE=7;

    co=[0.8500    0.3250    0.0980;         
        0         0.4470    0.7410;    
        0.9290    0.6940    0.1250;
        0.4940    0.1840    0.5560;
        0.4660    0.6740    0.1880;
        0.3010    0.7450    0.9330;
        0.6350    0.0780    0.1840
        0.4470    0.7410    0
        0.7410    0         0.4470
        0         0         0];
    markers = {'o','s','d','^','>','v','<','p','h','x'};
    t=1;
    pause(1);figure;pause(1)
    set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])

    for e=1:E
        if strcmp(expN,'05');lab{e}=cell2mat(recTyp(e,1));
        elseif strcmp(expN,'06');lab{e}=strcat(cell2mat(recTyp(e,1)),sprintf(' ($L=%s$)',cell2mat(recTyp(e,2)))); 
        elseif strcmp(expN,'07');lab{e}=strcat(cell2mat(recTyp(e,1)),sprintf(' ($L=%s$)',cell2mat(recTyp(e,2))));
        end
    end

    %Data ranges
    rangY=[inf 0];
    rangX=[0 nExtern];
    rangY=repmat(rangY,[S 1]);
    for s=1:S
        for v=1:V        
            for e=1:E
                rangY(s,1)=min(rangY(s,1),min(errF{s}{v+dfv}{e}));
                rangY(s,1)=min(rangY(s,1),min(errFGTNOTR{s}{v+dfv}{e}(:,1)));
                rangY(s,2)=max(rangY(s,2),max(errF{s}{v+dfv}{e}));
                rangY(s,2)=max(rangY(s,2),max(errFGTNOTR{s}{v+dfv}{e}(:,1)));
            end
        end
    end

    if ~strcmp(expN,'07')
        rangYRef=rangY(:,1);
        for s=1:S
            rangY(s,1)=min(rangY(s,1),errFGTNO{s});
            rangYRef(s)=rangY(s,1);
            rangY(s,:)=rangY(s,:)/rangYRef(s);        
            rangY(s,1)=rangY(s,1)/1.5;
        end
    else
        rangYRef=rangY(:,1);
        rangY=bsxfun(@rdivide,rangY,rangYRef);
        rangY(:,1)=rangY(:,1)/1.5;
    end
    rangY(:,2)=max(rangY(:,2));
    rangY(:,1)=min(rangY(:,1));

    for s=1:S
        for v=1:V    
            subtightplot(S,V,v+(s-1)*V,[0.02 0.02],[0.07 0.05],[0.03 0.05])           
            limsaxis=[rangX rangY(s,:)];
            for e=1:E        
                semilogy(effIt{s}{v+dfv}{e},errF{s}{v+dfv}{e}/rangYRef(s),'Color',co(e,:),'Marker',markers{e},'LineWidth',2);
                hold on      
            end    
            for e=1:E
                for n=1:size(errFGTNOTR{s}{v+dfv}{e},1);semilogy(1:limsaxis(2),(errFGTNOTR{s}{v+dfv}{e}(n,1)/rangYRef(s))*ones(1,limsaxis(2)),'Color',co(e,:),'LineWidth',1,'LineStyle','--');end
            end
            if ~strcmp(expN,'07');semilogy(1:limsaxis(2),(errFGTNO{s}/rangYRef(s))*ones(1,limsaxis(2)),'Color',co(end,:),'LineWidth',1,'LineStyle','--');end
            axis(limsaxis)
            set(gca,'fontsize',FontSizeA);

            if s==S;xlabel('$j$','interpreter','latex','FontSize',FontSizeB);else set(gca,'XTickLabel',[]);end

            if v==1;ylabel('$r$','interpreter','latex','FontSize',FontSizeB,'Rotation',0);else set(gca,'YTickLabel',[]);end        

            %HERE IT MAY BE RANGY(1) FOR FIGN NOT EQUAL TO 8
            if strcmp(expN,'05');text(nExtern/10,rangY(s,1)*1.2,sprintf('$\\theta=$%d$^{\\circ}$ / $M=$%d',theta(1,v),NS(1,s)),'interpreter','latex','Fontsize',FontSizeB);
            elseif strcmp(expN,'06');text(nExtern/10,rangY(s,1)*1.2,sprintf('$\\theta=$%d$^{\\circ}$ / $M=$%d',theta(1,v),NS(1,s)),'interpreter','latex','Fontsize',FontSizeB);
            elseif strcmp(expN,'07');text(nExtern/10,rangY(s,1)*1.2,sprintf('$\\theta=$%d$^{\\circ}$ / $M=$%d / $R=$%.1f$\\times$%.1f',theta(1,v+1),NS(1,s),NS(2,s),NS(2,s)),'interpreter','latex','Fontsize',FontSizeB);
            end

            if s==1 && v==V
                AX=legend(lab,'Location','NorthEast','FontSize',FontSizeB,'interpreter','latex');    
                LEG = findobj(AX,'type','text');
                set(AX,'Position',get(AX,'Position')+[0.05 0.04 0 0])%To the right and to the top         
            end
            grid on
        end
    end

    if strcmp(expN,'05')
        titA='GT / Uncorrected Sequential / Corrected Sequential / Uncorrected Random-checkered / Corrected Random-checkered';
        titB={'GT','Uncorrected Sequential','Corrected Sequential','Uncorrected Random-checkered','Corrected Random-checkered'};
        %10deg, 64 shots (normal)       
        %5deg, 16 shots (quick)
        xFig{1}=xGTNO{2-quick}{1};
        xFig{2}=xNOTR{2-quick}{3-quick}{1};
        xFig{3}=xEst{2-quick}{3-quick}{1};
        xFig{4}=xNOTR{2-quick}{3-quick}{4};
        xFig{5}=xEst{2-quick}{3-quick}{4};
        pause(1);figure;pause(1)
        imshow(1000*abs(cat(2,xFig{:})),[0 25],'Border','tight')
        colorbar('FontSize',FontSizeD);
        set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])
        title(titA,'FontSize',FontSizeB-6)

        for n=1:5;xFig{n}=xFig{n}-xGT;end
        pause(1);figure;pause(1)
        imshow(1000*abs(cat(2,xFig{:})),[0 2.5],'Border','tight')
        colorbar('FontSize',FontSizeD);
        set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])
        title(titA,'FontSize',FontSizeB-6)

        for n=1:length(xFig);fprintf('SNR %s: %.2f dB\n',titB{n},10*log10(normm(xGT)/normm(xFig{n})));end
    end

    if strcmp(expN,'07')             
        if quick
            titA='R=1 No motion / R=1 Known motion / R=1 Estimated motion / R=1.6 No motion / R=1.6 Known motion / R=1.6 Estimated motion';
            titB={'R=1 No motion','R=1 Known motion','R=1 Estimated motion','R=1.6 No motion','R=1.6 Known motion','R=1.6 Estimated motion'};
        else
            titA='R=1 No motion / R=1 Known motion / R=1 Estimated motion / R=2 No motion / R=2 Known motion / R=2 Estimated motion';
            titB={'R=1 No motion','R=1 Known motion','R=1 Estimated motion','R=2 No motion','R=2 Known motion','R=2 Estimated motion'};
        end

        NL=length(xGTNOTR{1}{1}{1});
        %SENSE1/SENSE2-NoMotion/10deg-RandChecker (normal)
        %SENSE1/SENSE1.6-NoMotion/10deg-RandChecker (quick)
        xFig{1}=xGTNOTR{1}{1}{3}{NL};
        xFig{2}=xGTNOTR{1}{3}{3}{NL};
        xFig{3}=xEst{1}{3}{3};
        xFig{4}=xGTNOTR{2}{1}{3}{NL};
        xFig{5}=xGTNOTR{2}{3}{3}{NL};
        xFig{6}=xEst{2}{3}{3};
        pause(1);figure;pause(1)
        imshow(1000*abs(cat(2,xFig{:})),[0 25],'Border','tight')
        colorbar('FontSize',FontSizeD);
        set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])
        title(titA,'FontSize',FontSizeB-6)

        for n=1:length(xFig);xFig{n}=xFig{n}-xGT;end
        pause(1);figure;pause(1)
        imshow(1000*abs(cat(2,xFig{:})),[0 2.5],'Border','tight')
        colorbar('FontSize',FontSizeD);
        set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])
        title(titA,'FontSize',FontSizeB-6)

        for n=1:length(xFig);fprintf('SNR %s: %.2f dB\n',titB{n},10*log10(sum(abs(xGT(:)).^2)/sum(abs(xFig{n}(:)).^2)));end
    end    
    pause(1)    
end