function ROI=computeROI(x,ext,sy,ev,fu)

%COMPUTEROI computes a spatial ROI from a mask
%   ROI=COMPUTEROI(X,{EXT},{SY},{EV})
%   * X is the mask of the structure of interest
%   * {EXT} is the range of security of the ROI with respect to the mask, defaults to 0
%   * {SY} forces the generated ROI to be strictly symmetric with regard of the original
%   FOV of the data, i.e., same number of leading and trailing excluded samples, 
%   defaults to 0
%   * {EV} forces the resulting ROI size to be even. It partially overrides the SYM 
%   parameter so that the resulting ROI would come from an approximately symmetric 
%   sample exclusion, defaults to 0
%   * {FU} returns a N-D ROI potentially with N>3-D. Otherwise it returns a N<=3-D ROI
%   * ROI are the ranges of the computed ROI, different dimensions are arranged along the
%   rows while the columns describe (1) lower limit of the ROI, (2) upper limit of the 
%   ROI, (3) size of the original data, (4) number of elements in the ROI, (5) number of 
%   leading excluded elements, (6) number of trailing excluded elements
%

if ~exist('ext','var') || isempty(ext);ext=0;end
if ~exist('fu','var') || isempty(fu);fu=0;end
ND=numDims(x);
if ~fu;ND=min(ND,3);end
if ~exist('sy','var') || isempty(sy);sy=single(zeros(1,ND));end
if ~exist('ev','var') || isempty(ev);ev=single(zeros(1,ND));end

%INITIALIZE
x=single(abs(x)>1e-6);
N=size(x);
ROI=single(zeros(ND,6));
ROI(:,1)=1;ROI(:,2)=N(1:ND)';ROI(:,3)=N(1:ND)';

%DETECT ROI
a=[1 -1];
for l=1:2
    for n=1:ND
        while isempty(find(dynInd(x,ROI(n,l),n),1));ROI(n,l)=ROI(n,l)+a(l);end
    end
end

%FILL ROI ARRAY AND APPLY CONSTRAINS
ROI(:,1)=max(ROI(:,1)-ext,1);     
ROI(:,2)=min(ROI(:,2)+ext,ROI(:,3));  
ROI(:,4)=ROI(:,2)-ROI(:,1)+1;
ROI(:,5)=ROI(:,1)-1;
ROI(:,6)=ROI(:,3)-ROI(:,2);
for n=1:ND
    if sy(n)
        [~,indSm]=min(ROI(n,5:6));
        ROI(n,[3-indSm 4])=ROI(n,[3-indSm 4])+ROI(n,6)-ROI(n,5);
        ROI(n,6)=ROI(n,5);
    end
    if ev(n) && mod(ROI(n,4),2)~=0
        ROI(n,4)=ROI(n,4)-1;
        ROI(n,[1 5])=ROI(n,[1 5])+1;       
    end
end

