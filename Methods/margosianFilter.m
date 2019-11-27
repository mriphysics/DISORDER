function [x,indMB]=margosianFilter(x,Enc,e,di,filt,useh)

% MARGOSIANFILTER applies a Margosian-like filter to half scan acquisitions
%   X=MARGOSIANFILTER(X,ENC,{E},{DI},{FILT})
%   * X is the data to be filtered
%   * ENC is some encoding structure
%   * {E} is the echo/mix number
%   * {DI} is the filter direction: 1, the data is filtered (default) / 0, 
%   the filter is reversed
%   * {FILT} is the type of filter. Possibilities are 'ramp' and 'zefi'. It 
%   defaults to ramp
%   * {USEH} activates phase-constrain
%   ** X is the filtered data
%   ** INDMB contains the information about what indexes have been set to 0
%

if nargin<3 || isempty(e);e=1;end
if nargin<4 || isempty(di);di=1;end
if nargin<5 || isempty(filt);filt='ramp';end
if nargin<6 || isempty(useh);useh=0;end

gpu=isa(x,'gpuArray');
NDims=length(Enc.kRange);
N=size(x);N(end+1:NDims)=1;N(4:end)=[];

if Enc.FOVSize(3)==1 && N(3)>1;Enc.FOVSize(3)=N(3);end
kGridA=generateGrid(Enc.FOVSize,gpu,Enc.FOVSize,ceil((Enc.FOVSize+1)/2));
kGridB=generateGrid(N,gpu,N,ceil((N+1)/2));

invert=1e-3;
invert2=1e-6;if invert>=invert2;invert2=0;end
indMA=cell(NDims,3);indMB=cell(NDims,3);
for m=1:NDims
    NR=diff(Enc.kRange{m}(e,:))+1;
    if NR~=1 && NR~=Enc.AcqSize(m)
        pfFactor=NR/Enc.AcqSize(m);     
        cutOffInd=round(Enc.FOVSize(m)*(1-pfFactor)+1);        
        cutOff=gather(abs(kGridA{m}(cutOffInd)-0.5));        
        indMA{m}{1}=(kGridA{m}(:)<=-cutOff);indMA{m}{2}=(kGridA{m}(:)>-cutOff & kGridA{m}(:)<=cutOff);indMA{m}{3}=(kGridA{m}(:)>cutOff);
        regSamp=N(m)-Enc.FOVSize(m);
        assert(regSamp>=0,'Reconstruction spectrum is smaller than acquired spectrum');
        fillRegSamp=[ceil(regSamp/2) floor(regSamp/2)];
        if mod(Enc.FOVSize(m),2)==0;flip(regSamp);end
        indMB{m}{1}=vertcat(true(fillRegSamp(1),1),indMA{m}{1},true(fillRegSamp(2),1));
        indMB{m}{2}=vertcat(false(fillRegSamp(1),1),indMA{m}{2},false(fillRegSamp(2),1));
        indMB{m}{3}=vertcat(false(fillRegSamp(1),1),indMA{m}{3},false(fillRegSamp(2),1));
                
        kyRM=cell(1,3);for n=1:3;kyRM{n}=gather(kGridB{m}(indMB{m}{n}));end
        sz=ones(1,NDims);sz(m)=N(m);sz(end+1:2)=1;        
        H=single(zeros(sz));F=single(zeros(sz));
        
        if strcmp(filt,'ramp');H(indMB{m}{1})=invert;H(indMB{m}{2})=1+(1-invert)*kyRM{2}/cutOff;H(indMB{m}{3})=2-invert;
        else H(indMB{m}{1})=invert;H(indMB{m}{2})=1;H(indMB{m}{3})=1;end
        F(indMB{m}{2})=1-(abs(kyRM{2})/cutOff);
        
        if ~di;H=1./(H+invert2);end
        if gpu;H=gpuArray(H);F=gpuArray(F);end
        if abs(Enc.kRange{m}(e,1))>abs(Enc.kRange{m}(e,2))
            H=flip(H,m);F=flip(F,m);
            for n=1:3;indMB{m}{n}=flip(indMB{m}{n});end
        end        
        H=ifftshift(H,m);
        F=ifftshift(F,m);
        for n=1:3;indMB{m}{n}=ifftshift(indMB{m}{n});end
        y=filtering(x,F);
        x=filtering(x,H);
        if useh;x=real(x.*conj(sign(y)));end    
    end
end
