function [Tou,F,NE]=compressMotion(Tin,F,res,parXT)

%COMPRESSMOTION   Bins the motion states using a Haar wavelet 
%   decomposition for quicker reconstructions
%   [TOU,F,NE]=COMPRESSMOTION(TIN,F,RES,PARXT)
%   * TIN are the input transforms for the motion states
%   * F are the Fourier encoding matrices with the different motion
%   states arranged as cell column elements
%   * RES is the normalized, isotropic equivalent data resolution
%   * PARXT are the parameters for aligned reconstruction
%   ** TOU are the transforms for the compressed motion states
%   ** F are the Fourier encoding matrices with the different compressed
%   motion states arranged as cell column elements
%   ** NE are the starting indexes of the Fourier encoding matrices for the
%   different compressed motion states
%

%INITIALIZATION
NT=size(Tin);
NSt=length(F{1});%Number of states
Nre=NSt-NT(5);%Number of states not associated with motion but with boundary conditions of the spectrum (samples outside the elliptic shutter). Corresponds to the number of repeats

Tou=Tin;
%COMPRESSION
perm=[6 5 1 2 3 4];
if NT(5)>1
    %TRANSFORM COMPRESSION
    Test=permute(Tin,perm);
    NL=floor(log2(prod(NT(1:5))));%Number of levels of wavelet decomposition
    C=cell(1,NT(6));L=cell(1,NT(6));
    for m=1:NT(6)%Number of parameters (6)
        if m<4;th=parXT.traLimX(1)*res;else th=convertRotation(parXT.traLimX(2)*res,'deg','rad');end%Binning factor depends on resolution
        [C{m},L{m}]=wavedec(Test(m,:),NL,'haar');
        C{m}(abs(C{m})<th)=0;
        Test(m,:) = waverec(C{m},L{m},'haar');
    end
    Test=ipermute(Test,perm);    

    %SAMPLE GROUPING
    n=1;
    NB=0; 
    indGroup=single(ones([1 NT(5)]));
    while n<=NT(5)
        Tcur=dynInd(Test,n,5);
        m=1;
        while m+n<=NT(5)
            Tnex=dynInd(Test,n+m,5);
            if any(abs(Tnex(:)-Tcur(:))>1e-6);break;end
            m=m+1;
        end
        indGroup(n:n+m-1)=m;
        n=m+n;
        NB=NB+1;
    end
    
    %MOTION GROUPING    
    Tou=single(zeros([NT(1:4) NB NT(6)]));
    %Inner samples
    n=1;nb=1;     
    while n<=NT(5)
        nG=n:n+indGroup(n)-1;
        Tou=dynInd(Tou,nb,5,mean(dynInd(Test,nG,5),5));
        nGr=nG(nG>nb);        
        for f=1:2;F{f}{nb}=cat(1,F{f}{nG});F{f}(nGr)={[]};end
        n=n+indGroup(n);
        nb=nb+1;
    end    
    %Outer samples and clean non used motion states / harmonics
    nbG=nb:nb+Nre-1;
    nG=n:n+Nre-1;nGr=nG(nG>nb+Nre-1);
    for f=1:2
        F{f}(nbG)=F{f}(nG);F{f}(nGr)={[]};
        F{f}=F{f}(~cellfun(@isempty,F{f}));
    end
end

%INDEX GROUPING
NE=single(zeros(1,length(F{1})+1));
for n=1:length(F{1});NE(n+1)=size(F{1}{n},1);end
NE=cumsum(NE);
