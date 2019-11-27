function visResiduals(We,outlWe,t,Residuals,ktraj,pau,folderName,fileName)

%VISRESIDUALS   Visualizes reconstruction residuals
%   VISRESIDUALS({WE},{OUTLWE},{T},{RESIDUALS},{KTRAJ},PAU,{FOLDERNAME},{FILENAME})
%   * {WE} is the inverse of the median error per motion state normalized 
%   to the median among states
%   * {OUTLWE} are the detected outliers in the motion states
%   * {T} is the time of the motion states
%   * {RESIDUALS} are the residuals in the spectrum
%   * {KTRAJ} are the trajectories
%   * {PAU} indicates whether to pause the execution, it defaults to 1
%   * {FOLDERNAME} gives a folder where to write the results
%   * {FILENAME} gives a file where to write the results
%

if nargin<1;We=[];end
if nargin<2;outlWe=[];end
if nargin<3;t=[];end
if nargin<4;Residuals=[];end
if nargin<5;ktraj=[];end
if nargin<6 || isempty(pau);pau=1;end
if nargin<7;folderName=[];end
if nargin<8;fileName=[];end

if ~iscell(folderName);folderName{1}=folderName;folderName{2}=folderName{1};end

co=[     0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330];

FontSizeA=30;
FontSizeB=24;
LineWidth=3;
MarkerSize=20;

if ~isempty(We)
    figure
    if ~isempty(t);xlab='$t$ (s)';else xlab='$t$ (au)';end
    if isempty(t);t=1:NS;end
    plot(t,We(:),'Color',co(1,:),'LineWidth',LineWidth,'LineStyle','-')
    hold on
    plot(t(outlWe),We(outlWe),'*','Color',co(2,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    if t(end)>t(1);xlim([t(1) t(end)]);end
    xlabel(xlab,'Interpreter','latex','Color',[1 1 1],'FontSize',FontSizeA+6)
    grid on        
    set(gca,'Color','none','XColor',[1 1 1],'YColor',[1 1 1],'FontSize',FontSizeA)
    dire{1}='Weights';
    if any(outlWe(:)~=0);dire{2}='Outliers';end
    AX=legend(dire);
    LEG = findobj(AX);
    set(LEG,'Color','none','TextColor',[1 1 1],'Location','NorthWest','FontSize',FontSizeA)          
    set(gcf,'Color',[0 0 0])
    set(gcf, 'Position', get(0,'Screensize'))  
    if pau==2 && ~isempty(folderName{1}) && ~isempty(fileName)%WE SIMPLY WRITE TO FILE
        if ~exist(folderName{1},'dir');mkdir(folderName{1});end
        %print(strcat(folderName,filesep,fileName),'-dpng');
        export_fig(strcat(folderName{1},filesep,fileName,'.png'));
        close all
    end
end
if ~isempty(Residuals)
    figure
    if ~isempty(ktraj)
        kMin=min(ktraj,[],1);
        kMax=max(ktraj,[],1);        
        for n=1:2;kGrid{n}=kMin(n):kMax(n);end
        kGrid{1}=kGrid{1}';
        kGrid{2}=repmat(kGrid{2},[1 size(Residuals,3)]);    
        for m=1:2;Residuals=fftshift(Residuals,m);end
        Residuals=log(Residuals(:,:));
        imagesc('XData',kGrid{2},'YData',kGrid{1},'CData',Residuals);
    else
        imagesc('CData',Residuals);
    end
    colormap(jet)
    axis image
    title('log($r$)','Interpreter','latex','FontSize',FontSizeB)
    xlabel('$k_3$','Interpreter','latex','FontSize',FontSizeB)
    ylabel('$k_2$','Interpreter','latex','FontSize',FontSizeB,'Rotation',0)
    if pau==2 && ~isempty(folderName{2}) && ~isempty(fileName)%WE SIMPLY WRITE TO FILE
        if ~exist(folderName{2},'dir');mkdir(folderName{2});end
        %print(strcat(folderName,filesep,fileName),'-dpng');
        export_fig(strcat(folderName{2},filesep,fileName,'.png'));
        close all
    end
end
if pau==1;pause;end
