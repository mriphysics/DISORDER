function visMotion(rec,voxSiz,t,pau,folderName,fileName,yl,xl,op)

%VISMOTION   Visualizes motion estimates
%   VISMOTION(REC,{VOXSIZ},{PAU},{FOLDERNAME},{FILENAME})
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   * {VOXSIZ} is the spacing of the image space
%   * {T} is the time of the motion states
%   * {PAU} serves to pause the execution, it defaults to 1
%   * {FOLDERNAME} gives a folder where to write the results
%   * {FILENAME} gives a file where to write the results
%   * {YL} are the limits for translations (first row) and rotations 
%   (second row)
%   * {XL} are the temporal limits
%   * {OP} serves to set a different opacity for different segments in the
%   xaxis
%

if nargin<3;t=[];end
if nargin<4 || isempty(pau);pau=1;end
if nargin<5;folderName=[];end
if nargin<6;fileName=[];end
if nargin<7;yl=[];end
if nargin<8;xl=[];end
if nargin<9;op=[];end

co=[     0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330];

FontSizeA=30;
ModLab=15;
LineWidth=6;

if pau==2;invCol=1;else invCol=0;end

if ~isstruct(rec);T=rec;rec=[];rec.T=T;end%In case that the transforms are input directly

if ~isempty(t);xlab='$t$ (s)';else xlab='$t$ (au)';end

if isfield(rec,'T')
    T=rec.T;
    ND=ndims(T);NT=size(T);
    NS=prod(NT(1:ND-1));
    T=reshape(T,[NS NT(ND)]);
    if isempty(t);t=1:NS;end
    if exist('voxSiz','var') && ~isempty(voxSiz);T(:,1:3)=bsxfun(@times,T(:,1:3),voxSiz);un='mm';else un='pix';end
    T(:,4:6)=180*T(:,4:6)/pi;
    figure

    %if ~isempty(xl)
    %    if size(t,3)==2;t=cat(3,dynInd(t,1,3),mean(t,3),dynInd(t,2,3));end
    %end
    if ~isempty(op) && size(t,3)==3;t=cat(3,dynInd(t,1,3),dynInd(t,3,3));end

    if ~isempty(op) && size(t,3)==2
        for l=1:2*size(t,1)
            if l<=size(t,1)+1;mark='-';else mark=':';end
            if l<=size(t,1)+1;LineW=LineWidth;else LineW=2*LineWidth/3;end
            if l>1;yyaxis left;end
            for n=1:3
                if l==1
                    tt=nan*ones(1,2);
                    TT=nan*ones(1,2);
                    OO=1;
                elseif l<=size(t,1)+1
                    tt=t(l-1,:,:);tt=tt(:)';
                    TT=repelem(T(l-1,n),2);TT=TT(:)';
                    OO=op(l-1);               
                else
                    tt=cat(2,t(l-size(t,1)-1,1,2),t(l-size(t,1),1,1));
                    TT=cat(2,T(l-size(t,1)-1,n),T(l-size(t,1),n));
                    OO=min(op(l-size(t,1)-1:l-size(t,1)));
                end         
                plot(tt,TT,mark,'Color',[co(n,:) OO],'LineWidth',LineW);%,'LineStyle','-')
                hold on
            end            

            if l==2
               set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSizeA)
            end

            
            if l==2*size(t,1)
                if ~isempty(xl);xlim(xl);else xlim([t(1) t(end)]);end
                if ~isempty(yl);ylim(yl(1,:));end
                xlabel(xlab,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
                ylabel(sprintf('Translation (%s)',un),'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
                grid on
            end
           
            if l>1;yyaxis right;end
            for n=4:6
                if l==1
                    tt=nan*ones(1,2);
                    TT=nan*ones(1,2);
                    OO=1;
                elseif l<=size(t,1)+1
                    tt=t(l-1,:,:);tt=tt(:)';
                    TT=repelem(T(l-1,n),2);TT=TT(:)';
                    OO=op(l-1);
                else
                    tt=cat(2,t(l-size(t,1)-1,1,2),t(l-size(t,1),1,1));
                    TT=cat(2,T(l-size(t,1)-1,n),T(l-size(t,1),n));
                    OO=min(op(l-size(t,1)-1:l-size(t,1)));
                end     
                plot(tt,TT,mark,'Color',[co(n,:) OO],'LineWidth',LineW);%,'LineStyle','-')
                hold on
            end

            if l==2
               set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSizeA)
            end

            if l==2*size(t,1)
                if ~isempty(xl);xlim(xl);else xlim([t(1) t(end)]);end    
                if ~isempty(yl);ylim(yl(2,:));end
                xlabel(xlab,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
                ylabel('Rotation (deg)','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
            end
        end
    else
        yyaxis left
        set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSizeA)
        
        yyaxis right
        set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSizeA)

        yyaxis left               
        for n=1:3
            if size(t,3)==2                        
                tt=t;
                tt=repmat(tt,[1 2 1]);
                tt=permute(tt,[2 3 1]);
                tt=tt(:);
                TT=repelem(T(:,n),4);
                TT(1:4:end)=0;
                TT(4:4:end)=0;
                plot(tt,TT,'Color',co(n,:),'LineWidth',LineWidth,'LineStyle','-')
            elseif size(t,3)==3
                %tt=t;
                %dt=tt(min(2,size(tt,1)),1,3)-tt(1,1,3);
                %tt=tt(1,1,1)+dt*(0:size(tt,1)-1);
                %tt=repelem(tt,2);
                %TT=repelem(T(:,n),2);
                %TT=[TT(1);TT(1:end-1)];
                tt=t(:,:)';tt=tt(:);       
                TT=repelem(T(:,n),3)';
                plot(tt,TT,'Color',co(n,:),'LineWidth',LineWidth,'LineStyle','-')
            else
                plot(t,T(:,n),'Color',co(n,:),'LineWidth',LineWidth,'LineStyle','-')
            end
            hold on
        end
        if ~isempty(xl);xlim(xl);else xlim([t(1) t(end)]);end
        if ~isempty(yl);ylim(yl(1,:));end
        xlabel(xlab,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
        ylabel(sprintf('Translation (%s)',un),'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
        grid on   
        
        yyaxis right
        for n=4:6
            if size(t,3)==2
                tt=t;
                tt=repmat(tt,[1 2 1]);
                tt=permute(tt,[2 3 1]);
                tt=tt(:);
                TT=repelem(T(:,n),4);
                TT(1:4:end)=0;
                TT(4:4:end)=0;
                plot(tt,TT,'Color',co(n,:),'LineWidth',LineWidth,'LineStyle','-')
            elseif size(t,3)==3
                %tt=t;
                %dt=tt(min(2,size(tt,1)),1,3)-tt(1,1,3);
                %tt=tt(1,1,1)+dt*(0:size(tt,1)-1);
                %tt=repelem(tt,2);
                %TT=repelem(T(:,n),2);
                %TT=[TT(1);TT(1:end-1)];
                tt=t(:,:)';tt=tt(:);       
                TT=repelem(T(:,n),3);
                plot(tt,TT,'Color',co(n,:),'LineWidth',LineWidth,'LineStyle','-')
            else                        
                plot(t,T(:,n),'Color',co(n,:),'LineWidth',LineWidth,'LineStyle','-')
            end
            hold on
        end
        if ~isempty(xl);xlim(xl);else xlim([t(1) t(end)]);end    
        if ~isempty(yl);ylim(yl(2,:));end
        xlabel(xlab,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
        ylabel('Rotation (deg)','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
    end
    grid on           
    
    if isfield(rec,'Par')
        if strcmp(rec.Par.Labels.FatShiftDir,'F') || strcmp(rec.Par.Labels.FatShiftDir,'H');dire{3}='FH';elseif strcmp(rec.Par.Labels.FatShiftDir,'A') || strcmp(rec.Par.Labels.FatShiftDir,'P');dire{3}='AP';else dire{3}='LR';end
        if strcmp(rec.Par.Labels.FoldOverDir,'HF') || strcmp(rec.Par.Labels.FoldOverDir,'FH');dire{2}='FH';elseif strcmp(rec.Par.Labels.FoldOverDir,'PA') || strcmp(rec.Par.Labels.FoldOverDir,'AP');dire{2}='AP';else dire{2}='LR';end
        if (strcmp(dire{3},'FH') || strcmp(dire{3},'AP')) && (strcmp(dire{2},'FH') || strcmp(dire{2},'AP'));dire{1}='LR';end
        if (strcmp(dire{3},'FH') || strcmp(dire{3},'LR')) && (strcmp(dire{2},'FH') || strcmp(dire{2},'LR'));dire{1}='AP';end
        if (strcmp(dire{3},'LR') || strcmp(dire{3},'AP')) && (strcmp(dire{2},'LR') || strcmp(dire{2},'AP'));dire{1}='FH';end
        for m=1:3;dire{m+3}=dire{m};end
        aux=dire{4};dire{4}=dire{6};dire{6}=aux;   
        for m=1:3;dire{m}=strcat('Tra-',dire{m});end
        for m=4:6;dire{m}=strcat('Rot-',dire{m});end
        AX=legend(dire);
        LEG = findobj(AX);
        set(LEG,'Color','none','TextColor',[1 1 1]*(1-invCol),'Location','NorthWest','FontSize',FontSizeA)            
    end  
    
    set(gcf,'Color',[0 0 0]+invCol)
    
    set(gcf, 'Position', get(0,'Screensize'))     
end
if pau==1;pause;end
if pau==2 && ~isempty(folderName) && ~isempty(fileName)%WE SIMPLY WRITE TO FILE
    if ~exist(folderName,'dir');mkdir(folderName);end
    %print(strcat(folderName,filesep,fileName),'-dpng');
    export_fig(strcat(folderName,filesep,fileName,'.png'));
    close all
end
