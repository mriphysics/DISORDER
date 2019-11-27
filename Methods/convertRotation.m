function x=convertRotation(x,inp,out)

%CONVERTROTATION   Converts units of rotation
%   X=CONVERTROTATION(X,INP,OUT)
%   * X is a rotation in units given by INP
%   * INP are the input units of rotation (either 'rad' or 'deg')
%   * OUT are the output units of rotation (either 'rad' or 'deg')
%   ** X is a rotation in units given by OUT
%

if strcmp(inp,'rad') && strcmp(out,'deg');x=180*x/pi;
elseif strcmp(inp,'deg') && strcmp(out,'rad');x=pi*x/180;
end