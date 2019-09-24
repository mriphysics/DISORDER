function [ySt,FSt,GSt,iStInd,nSt,mSt,nSa]=buildHarmonicSampling(y,G,H,NStates)

%BUILDHARMONICSAMPLING generates the cell structures corresponding to a 
%given partition of the sampled space into motion states
%   [YST,FST,GST,IST,ISTMOD,NST,MST,NSA]=BUILDHARMONICSAMPLING(Y,G,H,NSTATES)
%   * Y is the sampled data
%   * G is a sampling filter
%   * H corresponds to the motion states
%   * NSTATES is the number of motion states
%   ** YST are the samples with the different motion states arranged as a
%   column vector
%   ** FST are the Fourier encoding matrices with the different motion
%   states arranged as cell column elements
%   ** GST are the filter weights with the different motion states 
%   arranged as cell column elements
%   ** IST are the spectral and repetition indexes with the different 
%   motion states arranged as cell elements
%   ** NST are the starting indexes of the Fourier encoding matrices for 
%   the different motion states
%   ** MST is an array with the motion states that each sample corresponds 
%   to
%   ** NSA are the number of samples for the different motion states
%

%INITIALIZATION
gpu=isa(y,'gpuArray');
NY=size(y);NY(end+1:5)=1;
NH=size(H);H(end+1:3)=1;
H=H(:);
F=buildStandardDFTM(NY(1:2),1,gpu);
FSt=cell(2,NStates+1);ySt=cell(1,NStates+1);GSt=cell(1,NStates+1);iStInd=cell(1,NStates+1);mSt=H(:);
y=permute(y,[1 2 5 3 4]);
y=reshape(y,[prod(NY([1:2 5])) 1 NY(3:4)]);
if ~isempty(G);G=G(:);end
c=0;
for s=1:NStates+1 
    iStInd{s}=find(H==mod(s,NStates+1));
    iStSub=ind2subV(NH,iStInd{s});    
    for m=1:2;FSt{m}{s}=F{m}(iStSub(:,m),:);end
    ySt{s}=dynInd(y,iStInd{s},1);
    mSt(c+1:c+length(iStInd{s}))=s;
    c=c+length(iStInd{s});
    if ~isempty(G);GSt{s}=G(iStInd{s});end
end
ySt=cat(1,ySt{:});
GSt=cat(1,GSt{:});
nSa=accumarray(mSt,ones(1,length(mSt))')';
nSa(end+1:NStates+1)=0;nSa=nSa(1:NStates);
nSt=[0 cumsum(nSa(1:end))];
nSt=[nSt nSt(end)+size(FSt{1}{NStates+1},1)];
