function [c,y]=solveC(x,y,E,EH)

%SOLVEC   Computes the chi2-factors for a given reconstruction according to
%[1] R Winkelmann, P Bornert, O Dossel, "Ghost artifact removal using a 
%parallel imaging approach," Magn Reson Med, 54:1002-1009, 2005.
%   Y=SOLVEC(X,Y,E,EH)
%   * X is the reconstructed data
%   * Y is the measured data
%   * E is an encoding structure
%   * EH is a decoding structure
%   ** C is the chi2-factor
%   ** Y are the residuals in k-space
%

y=encode(x,E)-y;
if isfield(E,'Sf');E=rmfield(E,'Sf');end
if isfield(E,'Tr');E=rmfield(E,'Tr');end
if isfield(E,'Zf');E=rmfield(E,'Zf');end
%if isfield(EH,'Mb');EH.Mb=sqrt(EH.Mb./(sqrt(multDimSum(EH.Mb.^2,4))/size(EH.Mb,4)));end
if isfield(EH,'Mb');EH=rmfield(EH,'Mb');end
c=decode(y,EH,E);
ND=numDims(c);
NC=size(c,4);
c=normm(c,[],4:max(ND,4))/NC;
