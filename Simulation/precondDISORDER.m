function [P,SH]=precondDISORDER(S,reg)

%PRECONDDISORDER   Builds the preconditioner
%   P=PRECONDDISORDER(S,{REG})  
%   * S are the sensitivities
%   * {REG} is a regularization for inversion, it defaults to 0.001
%   ** P is the preconditioner
%   ** SH are the conjugated sensitivity maps
%

if nargin<2;reg=0.001;end
P=1./(normm(S,[],4)+reg);
if nargout>=2;SH=conj(S);end