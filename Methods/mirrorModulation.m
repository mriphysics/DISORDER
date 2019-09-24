function x=mirrorModulation(x,dirs,c)

% MIRRORMODULATION mirror filter modulation
%   X=MIRRORMODULATION(X,DIRS,CENTER) modulates a filter so that the zero
%   frequency is shifted to the Nyquist frequency
%   * X is the array to be modulated
%   * {DIRS} are the directions to modulate, it defaults to all directions
%   * {C} specifies the origin of modulation as floor(size(x)/2)+1+center 
%   (default is all zeros)
%   ** X is the modulated array
%

isd=isa(x,'double');
gpu=isa(x,'gpuArray');
N=size(x);
if nargin<2;dirs=1:length(N);end
if nargin<3;c=zeros(1,length(N));end
c=floor(N/2)+1+c;
rGrid=generateGrid(N,gpu,N,c);
for s=1:length(dirs)
    m=(-1).^rGrid{dirs(s)};
    x=bsxfun(@times,x,m);
end
if isd;x=double(x);end

end