function x=precondition(x,A,qr,di)

%PRECONDITION   Preconditions the residuals of inversion
%   X=PRECONDITION(X,A,{QR},{DI})
%   * X are the residuals before preconditioning
%   * A is the structure for preconditioning
%   * {QR} indicates whether we are using the QR algorithm in which case we 
%   apply the square root of the preconditioner, it defaults to 0
%   * {DI} indicates the direction of preconditioning (defaults to forward,
%   1)
%   ** X are the residuals after preconditioning
%

if isempty(A);return;end
if nargin<3 || isempty(qr);qr=0;end
if nargin<4 || isempty(di);di=1;end


if isfield(A,'Se')
    if qr==1;A.Se=sqrt(A.Se);end
    if di==0;A.Se=1./A.Se;end
    x=bsxfun(@times,x,A.Se);
end
if isfield(A,'Ps')
    if qr==1;A.Ps=sqrt(A.Ps);end
    if di==0;A.Ps=1./A.Ps;end
    x=filtering(x,A.Ps,A.mi);
end

