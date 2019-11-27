function x=mirroring(x,mirror,di,ty)

% MIRRORING mirrors a given array so as to guarantee periodicity
%   X=MIRRORING(X,MIRROR,DI,{TY})
%   * X is the array to be mirrored
%   * MIRROR is a flag to select whether a given direction has to be 
%   mirrored
%   * DI is the direction of mirroring: 1->mirror / 0->demirror
%   * {TY} is the type of mirroring: 0->unsymmetric (2 replicas) / 
%   1->symmetric (3 replicas), it defaults to 0)
%   ** X is the mirrored image
%

if ~exist('ty','var') || isempty(ty);ty=0;end

N=size(x);
nDimsOu=length(mirror);
N(end+1:nDimsOu)=1;
nDimsIn=length(N);
assert(nDimsOu<=nDimsIn,'Mirroring dimensionality is larger than image dimensionality');
assert(ismember(ty,[0 1]),'Type of mirroring is %d and has to be either 0 or 1',ty);

for m=1:nDimsOu       
    if mirror(m)     
        if di
            v=1:N(m);vo=v;v=horzcat(v,flip(v));if ty;v=horzcat(flip(vo),v);end           
        else
            v=1+ty*N(m)/3:(1+ty)*N(m)/(2+ty);
        end
        x=dynInd(x,v,m);
    end
end

