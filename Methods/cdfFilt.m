function y=cdfFilt(x,ord,winSize,pad_opts)

% CDFFILT performs N order-statistic filtering on ND data. It is a ND 
% extension of https://uk.mathworks.com/matlabcentral/fileexchange/22044-ordfilt3
%   Y=CDFFILT(X,{ORD},{WINSIZE},{PAD_OPTS}) 
%   * X is the input data
%   * {ORD} is the filtering type. Options are: 'max': maximum, 'min':
%   minimum, 'med': median, 'mea': mean, 'sum': sum, 'std': standard
%   deviation, 'mad': mean absolute deviation, 'iqr': interquantile range,
%   'var': variance. If a number o, it returns the o-th element when 
%   sorting. It defaults to 'med'
%   * {WINSIZE} is the window diameter, it defaults to 3 in all dimensions
%   * {PAD_OPTS} is the used padding. Options are those accepted by the
%   padarray function. It defaults to 'replicate'
%   ** Y is the filtered data
%

ND=numDims(x);
if nargin<2 || isempty(ord);ord='med';end
if nargin<3 || isempty(winSize);winSize=3*ones(1,ND);end
if nargin<4 || isempty(pad_opts);pad_opts='replicate';end
if numel(winSize)==1;winSize=winSize*ones(1,ND);else winSize(end+1:ND)=1;end
assert(all(mod(winSize,2)==1),'Window size must be odd');

if isnumeric(ord)%Parse ord argument
    assert(ord>=1 && ord<=prod(winSize),'Order parameter (%d) is out of bounds ([1 %d])',ord,winSize^ND);
    if ord==1;ord='min';end%more efficient to use 'min' rather than sorting
    if ord==prod(winSize);ord='max';end%more efficient to use 'max' rather than sorting
end
assert(~ischar(ord) || strcmp(ord,'max') || strcmp(ord,'min') || strcmp(ord,'med') || strcmp(ord,'mea') || strcmp(ord,'sum') || strcmp(ord,'std') || strcmp(ord,'mad') || strcmp(ord,'iqr') || strcmp(ord,'var'),'Invalid ord parameter %s',ord);%check for valid ord string

w=(winSize-1)/2;%Radious
x=padarray(x,w,pad_opts);%Padded volume
M=1e7;%Assumed amount of free contiguous memory (datatype not considered!)
y=ordfilt_compute(x,ord,M,w,ND);%Do the recursive call:

function D=ordfilt_compute(x,ord,M,w,ND)%ordfilt_compute: Recursively called for computing the ND order-statistics:

% ORDFILT_COMPUTE is recursive called for computing the ND
% order-statistics
%   D=ORDFILT_COMPUTE(X,ORD,M,W,ND) 
%   * X is the input data
%   * ORD is the filtering type. Options are: 'max': maximum, 'min':
%   minimum, 'med': median, 'mea': mean, 'sum': sum, 'std': standard
%   deviation, 'mad': mean absolute deviation, 'iqr': interquantile range,
%   'var': variance. If a number o, it returns the o-th element when 
%   sorting
%   * M is the assumed amount of freee contiguous memory (datatype not
%   considered!)
%   * W is the window radious
%   * ND is the number of dimensions
%   ** D is the filtered data
%

S=size(x);S(end+1:ND)=1;
S=S(1:ND);
mid=round((S)/2);
if prod(2*w+1)*prod(S)>M%test if we should evaluate data block stats: 
    %No: split the volume (2^N sub-volumes and call ordfilt_compute on each)
    v=cell(1,ND);for n=1:ND;v{n}={1:mid(n)+w(n),mid(n)-w(n)+1:S(n)};end     

    N=2*ones(1,ND);
    ind=1:prod(N);
    sub=ind2subV(N,ind);

    A=cell(ND,2);
    for n=ind
        u=cell(1,ND);for l=1:ND;u{l}=v{l}{sub(n,l)};end
        z=dynInd(x,u,1:ND);
        assert(any(size(z)~=size(x)),'Not enough free contiguous memory to compute order statistics');%Not enough free memory, since block size has not been reduced
        A{1}{sub(n,1)}=ordfilt_compute(z,ord,M,w,ND);
        l=1;
        while sub(n,l)==2 && l<ND            
            if l<ND
                A{l+1}{sub(n,l+1)}=cat(l,A{l}{:});
                l=l+1;            
            end
        end
    end
    D=cat(ND,A{ND}{:});                       
else
    %Yes: compute it concatinated neighbouring elements along ND+1th dimension:    
    N=(2*w+1);
    Bfull=zeros([S-2*w prod(N)],'like',x); 
    ind=1:prod(N);
    sub=bsxfun(@minus,ind2subV(N,ind),w(:)')-1;
    for cc=ind
        indI=cell(1,ND);for l=1:ND;indI{l}=(w(l)+1:S(l)-w(l))+sub(cc,l);end                        
        Bfull=dynInd(Bfull,cc,ND+1,dynInd(x,indI,1:ND));
    end
    if isnumeric(ord);D=dynInd(sort(Bfull,ND+1),ord,ND+1);
    elseif strcmp(ord,'max');D=max(Bfull,[],ND+1);
    elseif strcmp(ord,'min');D=min(Bfull,[],ND+1);
    elseif strcmp(ord,'med');D=median(Bfull,ND+1);
    elseif strcmp(ord,'mea');D=mean(Bfull,ND+1);
    elseif strcmp(ord,'sum');D=sum(Bfull,ND+1);
    elseif strcmp(ord,'std');D=std(Bfull,0,ND+1);
    elseif strcmp(ord,'mad');D=median(abs(bsxfun(@minus,Bfull,median(Bfull,ND+1))),ND+1);%D=mad(Bfull,1,ND+1);
    elseif strcmp(ord,'iqr');D=iqr(Bfull,ND+1);
    elseif strcmp(ord,'var');D=var(Bfull,0,ND+1);
    end
end
