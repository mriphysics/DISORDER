function x=constrain(x,C)

%CONSTRAIN   Constrains a reconstruction solution according to C
%   X=CONSTRAIN(X,C)
%   * X is the reconstruction before applying the constrains
%   * C is the structure of constrains
%   ** X is the reconstruction solution after applying the constrains
%

if isempty(C);return;end

if isfield(C,'Ma');x=bsxfun(@times,x,C.Ma);end
if isfield(C,'Fi');x(C.Fi.Ma)=C.Fi.Va(C.Fi.Ma);end
if isfield(C,'Re') && C.Re;x=real(x);end
if isfield(C,'Po') && C.Po;x=max(x,0);end
if isfield(C,'mV');x=bsxfun(@minus,x,multDimMea(x,C.mV));end
if isfield(C,'mD');x=bsxfun(@minus,x,multDimMed(x,C.mD));end
