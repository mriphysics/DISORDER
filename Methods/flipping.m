function x=flipping(x,H)

% FLIPPING flips an image using a given center
%   X=FLIPPING(X,H)
%   * X is the array to be flipped
%   * H is a cell array with the centers along different dimensions, it
%   corresponds to index-one centers, i.e. in pixel coordinates
%   ** X is the flipped image
%

ND=length(H);
for n=1:ND
    if ~isempty(H{n});H{n}=-H{n}+0.5;end
end
x=shifting(x,H);
for n=1:ND
    if ~isempty(H{n})
        H{n}=-H{n};
        x=flip(x,n);
    end
end
x=shifting(x,H);