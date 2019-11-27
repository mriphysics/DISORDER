function fig0304

%FIG0304    Generates Figs. 3 and 4 of the manuscript "Motion corrected 
%magnetic resonance imaging with DISORDER: Distributed and Incoherent 
%Sample Orders for Reconstruction Deblurring using Encoding Redundancy, 
%L Cordero-Grande, G Ferrazzi, RPAG Teixeira, AN Price, J O’Muircheartaigh, 
%and JV Hajnal, arXiv:1910.00540, 2019.
%   FIG0304
%

pathData='../DISORDERData';%Path with data, modify if required

addpath(genpath('..'));

for p=1:2
    FontSizeA=10;
    FontSizeB=18;
    FontSizeC=28;

    S=30;
    E=121;
    kInit=zeros(S*E,2);
    fid=fopen(sprintf('%s/exampleSpectrum.txt',pathData),'r');%This is an exemplary manufacturer sampling order, roughly corresponding to the ZigZag Sequential case below, as it comes from a shot-based sequence
    for n=1:S*E;kInit(n,:)=fscanf(fid,'%d',2);end
    fclose(fid);
    fprintf('Number of echoes: %d\n',E);fprintf('Number of shots: %d\n',S);

    U=[5 6];
    if p==1
        SampleOrderV={'ZigZag'};
        SegmentOrderV={'Sequential','Checkered','RandomCheckered'};
    else
        SampleOrderV={'AlternatingZigZag'};
        SegmentOrderV={'Sequential','Checkered','RandomCheckered','Random'};
    end
    t=1;
    for segmentOrder=1:length(SegmentOrderV)
        SegmentOrder=SegmentOrderV{segmentOrder};
        for sampleOrder=1:length(SampleOrderV)
            SampleOrder=SampleOrderV{sampleOrder};

            K=max(kInit,[],1)-min(kInit,[],1)+1;
            [k,NP,K,kGrid,ShotNumber,EchoNumber]=samplingDISORDER(kInit,U,SegmentOrder,SampleOrder,K);

            %Plot the results
            figure()  
            subtightplot(2,3,1,[0.16 0.04],[0.03 0.07],[0.05 0.05])
            imagesc('XData',kGrid{2},'YData',kGrid{1},'CData',EchoNumber)
            colormap(jet)
            colorbar        
            set(gca,'FontSize',FontSizeB)
            axis image
            xlabel('$k_2$','Interpreter','latex','FontSize',FontSizeC)
            ylabel('$k_3$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
            title('$e$','Interpreter','latex','FontSize',FontSizeC)

            subtightplot(2,3,4,[0.16 0.04],[0.07 0.03],[0.05 0.05])
            imagesc('XData',kGrid{2},'YData',kGrid{1},'CData',ShotNumber)
            colormap(jet)
            colorbar
            set(gca,'FontSize',FontSizeB)
            axis image
            xlabel('$k_2$','Interpreter','latex','FontSize',FontSizeC)
            ylabel('$k_3$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
            title('$s$','Interpreter','latex','FontSize',FontSizeC)    

            subtightplot(2,3,2,[0.12 0.06],[0.05 0.05],[0.05 0.05])
            plot(k(:,1))
            hold on
            set(gca,'FontSize',FontSizeB)
            xlabel('$p$','Interpreter','latex','FontSize',FontSizeC)
            ylabel('$k_3$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
            axis([1 NP kGrid{1}(1) kGrid{1}(end)])
            grid on

            subtightplot(2,3,5,[0.12 0.06],[0.08 0.02],[0.05 0.05])
            plot(k(:,2))
            hold on
            set(gca,'FontSize',FontSizeB)
            xlabel('$p$','Interpreter','latex','FontSize',FontSizeC)
            ylabel('$k_2$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
            axis([1 NP kGrid{2}(1) kGrid{2}(end)])
            grid on

            subtightplot(2,3,3,[0.12 0.05],[0.05 0.05],[0.05 0.05])
            plot(diff(k(:,1)))
            hold on
            set(gca,'FontSize',FontSizeB)
            xlabel('$p$','Interpreter','latex','FontSize',FontSizeC)
            ylabel('$\mbox{d}k_3$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
            axis([1 NP-1 -K(1)+1 K(1)-1])        
            grid on

            subtightplot(2,3,6,[0.12 0.05],[0.08 0.02],[0.05 0.05])    
            plot(diff(k(:,2)))
            hold on
            set(gca,'FontSize',FontSizeB)
            ylabel('$\mbox{d}k_2$','Interpreter','latex','FontSize',FontSizeC,'Rotation',0)
            xlabel('$p$','Interpreter','latex','FontSize',FontSizeC)
            axis([1 NP-1 -K(2)+1 K(2)-1])
            grid on
            set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
            h = annotation('textbox', [0.6 1 0 0],'String',sprintf('%s - %s',SegmentOrder,SampleOrder),'FontSize',FontSizeC,'FitBoxToText', 'on');
        end
    end
end