function [resPyr,L,estT,resIso]=pyramidPlan(voxSiz,resolMax,NWend,accel)

%PYRAMIDPLAN  Constructs a multirresolution plan on the basis of some 
%parameters
%   [RESPYR,L,ESTT,RESISO]=pyramidPlan(VOXSIZ,RESOLMAX,NWEND,ACCEL)
%   * VOXSIZ is the voxel size
%   * RESOLMAX is the starting resolution of the pyramid
%   * NWEND is the final number of within shot subdivisions
%   * ACCEL is the number of levels for which no motion estimates are
%   obtained
%   ** RESPYR is the spatial resolution of each level
%   ** L is the total number of subdivisions
%   ** ESTT indicates whether to estimate motion for each level
%   ** RESISO is the equivalent isotropic resolution normalized to 0.7, 
%   which is the maximum expected
%

Lr=round(log2(max(resolMax/min(voxSiz),1)))+1;%Number of spatial multiresolution levels
resPyr=1./flip(2.^(0:Lr-1));%Spatial multiresolution pyramid        
L=Lr+round(log2(NWend));%Spatio-temporal number of multiresolution levels
resPyr(Lr+1:L)=1;
estT=ones(1,L);%Flag to estimate motion 
estT(max(L-accel(1)+1,1):L)=0;%Not to estimate motion to accelerate
resIso=sqrt(2)*((prod(voxSiz).^(1/3))./resPyr);%Equivalent resolution of isotropic voxels, normalized to a resolution of 0.7, which is the maximum we expect

L=max(L-accel(2),1);
estT=estT(1:L);
resIso=resIso(1:L);
resPyr=resPyr(1:L);
