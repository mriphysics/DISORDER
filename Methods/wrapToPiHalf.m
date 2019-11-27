function x=wrapToPiHalf(x)

%WRAPTOPIHALF   Wraps angles to pi/2 
%   X=WRAPTOPIHALF(X)
%   * X is the angle to wrap
%   ** X is the wrapped angle
%

x=mod(x,pi);

ind=x>pi/2;
x(ind)=x(ind)-pi;