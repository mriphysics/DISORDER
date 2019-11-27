function [F,FH]=build1DFTM(N,fn,gpu,re)

%BUILD1DFTM   Builds a standard DFT matrix
%   [F,FH]=BUILD1DFTM(N,{FN},{GPU},{RE})
%   * N is the dimension of the space
%   * {FN} indicates whether to generate fully unitary Fourier matrices. It
%   defaults to 0
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) matrices (empty, default depending on machine)
%   * {RE} is a flag to indicate that they are to be applied to a real 
%   image (defaults to 0)
%   ** F is a discrete Fourier transform matrix
%   ** FH is a an inverse discrete Fourier transform matrix
%

if nargin<2 || isempty(fn);fn=0;end
if nargin<3 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
if nargin<4 || isempty(re);re=0;end


F=single(dftmtx(N));
if gpu;F=gpuArray(F);end
if fn;F=F/sqrt(N);end
if ~re && nargout>1
    FH=conj(F);
elseif re
    Nre=ceil((N+1)/2);
    F=F(1:Nre,:);
    if re==1
        %This gives a real-only Fourier representation- We need real DC-real 1st harmonic-imag 2st harmonic-...-real Nyquist component (only if even)       
        F=cat(3,real(F),imag(F));
        F=permute(F,[3 1 2]);
        F=reshape(F,[],N);   
        if re==1
            F(2,:)=[];        
            if mod(N,2)==0;F(end,:)=[];end
        else
            F(2,:)=0;%Think this is not necessary
            if mod(N,2)==0;F(end,:)=0;end%Think this is not necessary
        end
        if nargout>1
            FH=2*F';    
            FH(:,1)=FH(:,1)/2;
            if mod(N,2)==0;FH(:,end-(re~=1):end)=FH(:,end-(re~=1):end)/2;end
        end    
    elseif re==2    
        FH=2*F';
        FH(:,1)=FH(:,1)/2;
        if mod(N,2)==0;FH(:,end)=FH(:,end)/2;end     
    end
end
if ~fn && nargout>1;FH=FH/N;end

