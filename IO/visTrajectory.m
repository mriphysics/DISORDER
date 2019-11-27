function visTrajectory(ktraj,pau,folderName,fileName)

%VISTRAJECTORY   Visualizes sampling trajectories
%   VISTRAJECTORY(KTRAJ,{PAU},{FOLDERNAME},{FILENAME})
%   * KTRAJ is the trajectory to be visualized
%   * {PAU} indicates whether to pause the execution, it defaults to 1
%   * {FOLDERNAME} gives a folder where to write the results
%   * {FILENAME} gives a file where to write the results
%

if nargin<2 || isempty(pau);pau=1;end
if nargin<3;folderName=[];end
if nargin<4;fileName=[];end

assert(numDims(ktraj)==3,'Number of dimensions for plotting trajectories should be 3 and it is %d',numDims(ktraj));
K=size(ktraj);
assert(K(3)==2,'Trajectories should be defined on a 2D space but they are on a %d space',K(3));

kMin=multDimMin(ktraj,1:2);
kMax=multDimMax(ktraj,1:2);
kSiz=kMax-kMin+1;kSiz=permute(kSiz,[1 3 2]);

EchoNumber=single(zeros(kSiz));
ShotNumber=EchoNumber;
NP=prod(K(1:2));

for e=1:K(1)
    for s=1:K(2)
        EchoNumber(ktraj(e,s,1)-kMin(1)+1,ktraj(e,s,2)-kMin(2)+1)=e;
        ShotNumber(ktraj(e,s,1)-kMin(1)+1,ktraj(e,s,2)-kMin(2)+1)=s;
    end
end
ktraj=reshape(ktraj,[prod(K(1:2)) K(3)]);

for n=1:2;kGrid{n}=kMin(n):kMax(n);end
kGrid{1}=kGrid{1}';

FontSizeA=10;
FontSizeB=18;
FontSizeC=28;

%To plot the results
figure

subtightplot(2,3,1,[0.16 0.04],[0.03 0.07],[0.05 0.05])
imagesc('XData',kGrid{2},'YData',kGrid{1},'CData',EchoNumber)
colormap(jet)
colorbar        
set(gca,'FontSize',FontSizeB)
axis image
xlabel('$k_3$','Interpreter','latex','FontSize',FontSizeC)
ylabel('$k_2$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
title('$e$','Interpreter','latex','FontSize',FontSizeC)

subtightplot(2,3,4,[0.16 0.04],[0.07 0.03],[0.05 0.05])
imagesc('XData',kGrid{2},'YData',kGrid{1},'CData',ShotNumber)
colormap(jet)
colorbar
set(gca,'FontSize',FontSizeB)
axis image
xlabel('$k_3$','Interpreter','latex','FontSize',FontSizeC)
ylabel('$k_2$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
title('$s$','Interpreter','latex','FontSize',FontSizeC)    

subtightplot(2,3,2,[0.12 0.06],[0.05 0.05],[0.05 0.05])
plot(ktraj(:,1))
hold on
set(gca,'FontSize',FontSizeB)
xlabel('$p$','Interpreter','latex','FontSize',FontSizeC)
ylabel('$k_2$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
axis([1 NP kGrid{1}(1) kGrid{1}(end)])
grid on

subtightplot(2,3,5,[0.12 0.06],[0.08 0.02],[0.05 0.05])
plot(ktraj(:,2))
hold on
set(gca,'FontSize',FontSizeB)
xlabel('$p$','Interpreter','latex','FontSize',FontSizeC)
ylabel('$k_3$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
axis([1 NP kGrid{2}(1) kGrid{2}(end)])
grid on

K=max(ktraj,[],1)-min(ktraj,[],1)+1;

subtightplot(2,3,3,[0.12 0.05],[0.05 0.05],[0.05 0.05])
plot(diff(ktraj(:,1)))
hold on
set(gca,'FontSize',FontSizeB)
xlabel('$p$','Interpreter','latex','FontSize',FontSizeC)
ylabel('$\mbox{d}k_2$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
axis([1 NP-1 -K(1)+1 K(1)-1])        
grid on

subtightplot(2,3,6,[0.12 0.05],[0.08 0.02],[0.05 0.05])    
plot(diff(ktraj(:,2)))
hold on
set(gca,'FontSize',FontSizeB)
ylabel('$\mbox{d}k_3$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
xlabel('$p$','Interpreter','latex','FontSize',FontSizeC)
axis([1 NP-1 -K(2)+1 K(2)-1])
grid on
set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);


set(gcf, 'Position', get(0,'Screensize'))
if pau==1;pause;end
if pau==2 && ~isempty(folderName) && ~isempty(fileName)%WE SIMPLY WRITE TO FILE
    if ~exist(folderName,'dir');mkdir(folderName);end
    %print(strcat(folderName,filesep,fileName),'-dpng');
    export_fig(strcat(folderName,filesep,fileName,'_Trajectory.png'));
    close all
end
