function varargout=parUnaFun(x,f,varargin)

%PARUNAFUN   Calls a unary function for different datasets
%   VARARGOUT=PARUNAFUN(X,F,VARARGIN)
%   * X is a cell array containing the input datasets
%   * F is a function handle describing the operation to perform on those
%   datasets
%   * VARARGIN are the arguments of the operation to be performed
%   * VARARGOUT is a variable size output with number of variables given by
%   the number of datasets in the input cell array
%

NO=nargout;
NI=length(x);
assert(NO==NI,'Number of outputs (%d) not equal to number of inputs (%d)',NO,NI);

varargout=cell(1,NO);
for o=1:NO;varargout{o}=f(x{o},varargin{:});end