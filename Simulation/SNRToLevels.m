function l=SNRToLevels(l,xGT,S,di)

%SNRTOLEVELS   Computes the noise levels to be added to the samples for a
%prescribed average spatial SNR
%   L=SNRTOLEVELS(L,XGT,S,{DI}) computes noise levels from SNR or SNR
%   from noise levels. Note the SNR is defined as mean(abs(XGT).^2), noise
%   is mean(abs(XNO).^2) and noise is generated using the plugNoise 
%   function
%   * L is a vector with the input levels, in dB if DI=1 and in noise units
%   if di=0
%   * XGT is the ground truth image
%   * S is the coil-array sensitivity map
%   * {DI} is the direction of conversion, from dB to noise levels (1,
%   default) or from noise levels to dB, 0. NOTE: THE CASE DI=0 IS NOT
%   IMPLEMENTED YET!
%   ** L is a vector with the output levels, in dB if DI=0 and in noise
%   units if DI=1
%

if nargin<4 || isempty(di);di=1;end

%Average power
meanS=mean(abs(xGT(:)).^2);
%Average noise power for prescribed levels
meanN=meanS./(10.^(l/10));

%Reconstruction noise amplification
reg=1e-3;
P=1./(sum(abs(S).^2,4)+reg);
P=mean(P(:))/numel(xGT);

%Noise level
l=sqrt(meanN/(2*P));
