function [En,W]=computeEnergy(y,x,E,R,EH,no,di,non,de,avdi,sepnorm)

%COMPUTEENERGY   Computes the energy of the solution of a least squares 
%problem with quadratic regularization. It also serves to provide the
%weights for an IRWLS iteration
%   [EN,W]=COMPUTEENERGY(Y,X,E,{R},{EH},{NO},{DI},{NON},{DE},{AVDI},{SEPNORM})
%   * Y is the measured data
%   * X is the reconstructed data
%   * E is the encoding structure
%   * {R} is a regularization structure
%   * {EH} is the decoding structure (used for weighted least squares)
%   * {NO} is the norm of fitting to be optimized
%   * {DI} are the dimensions over which the previous norm applies (for 
%   other dimensions the norm is supposed to be 2)
%   * {NON} is the norm used for weight estimation (distinct from NO when
%   continuation is used)
%   * {DE} is the regularization for weight computation (its role comes 
%   from a Huber loss)
%   * {AVDI} avoids certain dimensions in the energy computation
%   * {SEPNORM} introduces a normalization separable along the dimensions
%   not included in the norm
%   ** EN is a row vector with the energy, with first component for the 
%   fidelity and second for the regularization
%   ** W is an array with the weights derived from the residuals
%

if nargin<6 || isempty(no);no=2;end
if nargin<7 || isempty(di);di=1:ndims(y);end
if nargin<8 || isempty(non);non=no;end
if nargin<9 || isempty(de);de=1e-9;end
if nargin<10;avdi=[];end
if nargin<11 || isempty(sepnorm);sepnorm=0;end

fi=y-encode(x,E);%Residuals in image space generally


if nargin>=5 && isfield(EH,'Mb');fi=bsxfun(@times,fi,sqrt(EH.Mb));end
if nargin>=5 && isfield(EH,'Mbe');fi=bsxfun(@times,fi,sqrt(EH.Mbe));end%Only to weight the energy but not to be applied when decoding

fi=abs(fi);
if ~isempty(avdi);En=multDimSum(fi.^no,setdiff(di,avdi));W=[];return;end

if no>-1
    if no~=2       
        W=max(de,fi).^((non-2)/2);
        if ~sepnorm
            V=(multDimSum(max(de,fi).^non,di)).^((2-non)/2);        
            W=bsxfun(@times,W,V);%"Normalization"
        end
        W=W.^2;
        if sepnorm;W=bsxfun(@rdivide,W,multDimSum(W,di));end        
    end
    
    fi=multDimSum(fi.^no,di);
    fi=fi.^(1/no);
    fi=gather(normm(fi));
else%HUBER LOSS    
    sigma=1.4826*median(fi(:));
    tau=1.345*sigma;
    fi2=fi<=tau;
    fi1=fi>tau;
    
    p=-1-no;%Huber, when 0 (no=-1), linear Huber when 1 (no=-2), Huber with total saturation    
    W=zeros(size(fi),'like',fi);
    W(fi2)=1;
    W(fi1)=1*(1-p)*tau^(1+p)./(fi(fi1).^(1+p));
    
    fi=gather(sum(fi(fi2).^2)+sum(2*(tau^(1+p))*(fi(fi1).^(1-p))-tau*tau));
end    

if nargin>=4
    re=regularize(x,R,1);
    reaux=abs(sum(re(:).*conj(x(:))));
    re=regularize(x,R,4);
    re=gather(reaux+abs(sum(re(:))));
else
    re=0;
end
En=[fi;re];
