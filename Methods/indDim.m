function x=indDim(x,ind,dim,y,indb,indf)

%INDDIM   Indexes X along a given dimension DIM using the values in array
%IND which refer to indexes along that dimension for each of the array 
%elements at all or some of the other dimensions
%   X=INDDIM(X,IND,{DIM},{Y},{INDB},{INDF})
%   * X is the data array
%   * IND is the index array
%   * {DIM} is the dimension the index refers to (1 by default)
%   * {Y} is a second input from which all its elements are assigned to the
%   corresponding indexes in X
%   * {INDB} is the array of indexes before the given dimension (to
%   accelerate indexing)
%   * {INDF} is the array of indexes after the given dimension (to
%   accelerate indexing)
%   ** X is the output data array
%

if nargin<3 || isempty(dim);dim=1;end
gpu=isa(ind,'gpuArray');
ind=double(ind);

N=size(x);N(end+1:dim)=1;ND=length(N);
Nor=N(dim);
NI=size(ind);NI(end+1:ND)=1;
NI(dim)=1;N(dim)=1;
assert(length(NI)<=length(N),'The dimensions of the index array (%d) are bigger than the dimensions of the data array (%d)',length(NI),length(N));
%assert(all(ind(:)<=NP) && all(ind(:)>=1),'Indexes out of range');%Was taking too much time, commented
%assert(all(NI==N | NI==1),'Some index and data array dimensions not compatible');%Disabled for more flexible operation
ind=(ind-1)*prod(N(1:dim-1));
ind=repmat(ind,N./NI);
if nargin<5 || isempty(indb);indb=1:prod(N(1:dim-1));else indb=double(indb);end
if gpu;indb=gpuArray(indb);end
if ~isempty(N(1:dim-1));indb=reshape(indb,[N(1:dim-1) ones(1,ND-dim+1)]);end
if nargin<6 || isempty(indf);indf=0:prod(N(dim+1:end))-1;else indf=double(indf);end
if gpu;indf=gpuArray(indf);end
if ~isempty(N(dim+1:end));indf=reshape(indf,[ones(1,dim) N(dim+1:end)]);end
indf=indf*prod([N(1:dim-1) Nor]);
ind=bsxfun(@plus,ind,bsxfun(@plus,indb,indf));
if nargin>=4 && ~isempty(y);x(ind)=y;else x=x(ind);end
