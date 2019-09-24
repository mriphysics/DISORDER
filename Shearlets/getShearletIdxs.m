function idxs=getShearletIdxs(ND,K,full,complex,varargin)

%GETSHEARLETIDXS computes a index set describing a shearlet system
%   IDXS=GETSHEARLETIDXS(ND,K,{FULL},{VARARGIN}) computes the descriptive
%   indexes of the shearlet system
%   * ND are the number of dimensions of the system
%   * K are the levels of shearing occurring on each scale
%   * {FULL} determines whether a full shearlet system is computed or if
%   shearlets lying on borders of cones / pyramids are grouped. Defaults to
%   0
%   * {COMPLEX} indicates whether complex shearlets are generated, so both
%   sides of the spectrum are not treated symmetrically. Defaults to 0
%   * {VARARGIN} indicates possible restrictions by means of pairs:
%       - 'TypeRestriction1': Possible restrictions: 'cones', 'scales', 
%   'shearings' (2D) and 'pyramids', 'scales', 'shearings1', 'shearings2' 
%   (3D)
%       - ValueRestriction1: Numerical value or Array specifying a 
%   restriction. If the type of the restriction is 'scales' the value 1:2 
%   ensures that only indexes corresponding the shearlets on the first two
%   scales are computed
%   ** IDXS MX(ND+1) matrix, where each row describes one shearlet in the 
%   format [cone scale shearing] (2D) / [cone scale shearing1 shearing2] 
%   (3D)
%

if nargin<3 || isempty(full);full=0;end
if nargin<4 || isempty(complex);complex=0;end

idxs=[];
includeLP=1;

scales=1:length(K);
shearings=-2^(max(K)):2^(max(K));shearings1=shearings;shearings2=shearings;
cones=1:2;
pyramids=1:3;
if complex;halfs=1:2;else halfs=1;end

for v=1:2:length(varargin)
    includeLP = 0;
    if strcmp(varargin{v},'scales');scales=varargin{v+1};
    elseif strcmp(varargin{v},'shearings');shearings=varargin{v+1}; 
    elseif strcmp(varargin{v},'shearings1');shearings1=varargin{v+1}; 
    elseif strcmp(varargin{v},'shearings2');shearings2=varargin{v+1}; 
    elseif strcmp(varargin{v},'cones');cones=varargin{v+1}; 
    elseif strcmp(varargin{v},'pyramids');pyramids=varargin{v+1}; 
    elseif strcmp(varargin{v},'halfs');halfs=varargin{v+1}; 
    end
end

if ND==2
    for cone=intersect(1:2,cones)
        for scale=intersect(1:length(K),scales)
            for half=intersect(1:2,halfs)
                for shearing=intersect(-2^K(scale):2^K(scale),shearings)
                    if full || cone == 1 || abs(shearing)<2^K(scale);idxs=vertcat(idxs,[cone scale half shearing]);end                        
                end
            end
        end
    end
    if includeLP || ismember(0,scales) || ismember(0,cones);idxs=vertcat(idxs,[0 0 0 0]);end
elseif ND==3
    for pyramid=intersect(1:3,pyramids)
        for scale=intersect(1:length(K),scales)
            for half=intersect(1:2,halfs)
                for shearing1=intersect(-2^K(scale):2^K(scale),shearings1)
                    for shearing2 = intersect(-2^K(scale):2^K(scale),shearings2)
                        if full || pyramid == 1 || (pyramid == 2 && abs(shearing2)<2^K(scale)) || (abs(shearing1)<2^K(scale) && abs(shearing2)<2^K(scale));idxs=vertcat(idxs,[pyramid scale half shearing1 shearing2]);end
                    end
                end
            end
        end
    end
    if includeLP || ismember(0,scales) || ismember(0,pyramids);idxs=vertcat(idxs,[0 0 0 0 0]);end
end
