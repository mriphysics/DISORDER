function visSegment(x,y,pau,orthv,sl,dy,str,folderName,fileName)

%VISSEGMENT   Visualizes a segmentation overlaid on a given image
%   VISSEGMENT(X,{Y},{PAU},{SL},{DY})
%   * X is the image
%   * {Y} is the segmentation
%   * {PAU} indicates whether to pause the execution, it defaults to 1
%   * {ORTHV} shows orthogonal views
%   * {SL} serves to select a given slice
%   * {DY} serves to select a given dynamic
%   * {STR} serves to write a title
%   * {FOLDERNAME} gives a folder where to write the results
%   * {FILENAME} gives a file where to write the results
%

N=size(x);N(end+1:5)=1;
if nargin<2;y=[];end
if nargin<3 || isempty(pau);pau=1;end
if nargin<4 || isempty(orthv);orthv=single(N(3)>1);end
if nargin<5 || isempty(sl);sl=mod(ceil(N(3)/2)+1,N(3))+1;end
if nargin<6 || isempty(dy);dy=1;end
if nargin<7;str='';end
if nargin<8;folderName=[];end
if nargin<9;fileName=[];end
assert(numel(dy)==1,'Dynamics are non-singleton');
assert(numel(sl)==1,'Slices are non-singleton');
N=[N(1:3) prod(N(4:end))];
x=reshape(x,N);
x=dynInd(x,dy,4);
if ~isempty(y)
    y=reshape(y,N);
    y=dynInd(y,dy,4);
end
if orthv
    [x,y]=extractOrthogonalPlanes(x,y);
else
    x=dynInd(x,sl,3);
    if ~isempty(y);y=dynInd(y,sl,3);end
end
x=double(gather(x));
x=abs(x);
x=x-min(x(:));
x=x/max(x(:));
if ~isempty(y)
    y=double(gather(y));
    y=imquantize(y,0.5:max(y(:)))-1;
end

figure
imshow(x,'Border','tight')
hold on
if ~isempty(y);contour(y,[0.5 0.5:max(y(:))],'EdgeColor',[0.4660 0.6740 0.1880],'LineWidth',2);end
set(gcf, 'Position', get(0,'Screensize'))
text(1,10,str,'FontSize',24,'Color',[0.8500 0.3250 0.0980],'Interpreter','latex');
if pau==1;pause;end
if pau==2 && ~isempty(folderName) && ~isempty(fileName)%WE SIMPLY WRITE TO FILE
    if ~exist(folderName,'dir');mkdir(folderName);end
    %print(strcat(folderName,filesep,fileName),'-dpng');
    export_fig(strcat(folderName,filesep,fileName,'.png'));
    close all
end
