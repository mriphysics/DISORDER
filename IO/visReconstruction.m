function visReconstruction(x,pau,sl,dy,xlab,ylab)

%VISRECONSTRUCTION   Visualizes a slice of each of the reconstructed 
%volumes
%   VISRECONSTRUCTION(X,{PAU},{SL},{DY},{XLAB},{YLAB})
%   * X is the reconstructed data
%   * {PAU} indicates whether to pause the execution, it defaults to 1
%   * {SL} serves to select a given slice
%   * {DY} serves to select a given dynamic
%   * {XLAB} serves to set the label of the x axis
%   * {YLAB} serves to set the label of the y axis
%

N=size(x);N(end+1:5)=1;
if nargin<2 || isempty(pau);pau=1;end
if nargin<3 || isempty(sl);sl=mod(ceil(N(3)/2)+1,N(3))+1;end
if nargin<4 || isempty(dy);dy=1;end
if nargin<5 || isempty(xlab);xlab='Phase encode';end
if nargin<6 || isempty(ylab);ylab='Readout';end
assert(numel(dy)==1,'Dynamics are non-singleton');
assert(numel(sl)==1,'Slices are non-singleton');
N=[N(1:3) prod(N(4:end))];
x=reshape(x,N);
x=dynInd(x,[sl dy],3:4);
FontSizeA=30;
figure
subtightplot(1,2,1,[0 0.075])
imshow(abs(x),[])
ylabel(ylab,'FontSize',FontSizeA);
xlabel(xlab,'FontSize',FontSizeA)
title('Magnitude','FontSize',FontSizeA)
subtightplot(1,2,2,[0 0.075])
imshow(angle(x),[-pi pi])
ylabel(ylab,'FontSize',FontSizeA);
xlabel(xlab,'FontSize',FontSizeA)    
title('Phase','FontSize',FontSizeA)
set(gcf, 'Position', get(0,'Screensize'))  
if pau;pause;end