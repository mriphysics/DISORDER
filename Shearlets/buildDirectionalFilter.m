function [h0,h1]=buildDirectionalFilter(DFtype,di,modu)

%BUILDDIRECTIONALFILTER builds a directional filter
%   [H0,H1]=BUILDDIRECTIONALFILTER(DFTYPE,{DI},{MODU}) builds a directional
%   filter based on the dfilters function of the Nonsubsampled Contourlet
%   Toolbox 
%   http://www.mathworks.de/matlabcentral/fileexchange/10049-nonsubsampled-contourlet-toolbox
%   To test the diamond filter pairs (for the PR condition for the FIR 
%   case), verify that conv2(h0,mirrorModulation(h1))+conv2(mirrorModulation(h0),h1)=2
%   (replace + with - for even size filters)
%   To test for orthogonal filter conv2(h,reverse2(h))+mirrorModulation(conv2(h,reverse2(h)))=2
%   * DFTYPE is the directional filter type
%   * {DI} indicates a decomposition (1) or reconstruction (0) filter. It 
%   defaults to  1
%   * {MODU} indicates whether to modulate the filter, it defaults to 1
%   ** H0 is the directional filter if MODU equals 1 or the lowpass of the
%   diamond filter pair
%   ** H1 is the high pass of the diamond filter pair in case MODU equals 0
%   and empty otherwise
%

if nargin<2 || isempty(di);di=1;end
if nargin<3 || isempty(modu);modu=1;end

t=[0 1 0;1 0 1;0 1 0]/4;% diamond kernel
switch DFtype
    case 'haar'
        h0=[1 1]/sqrt(2);
        h1=(-1)^di*[1 -1]/sqrt(2);	
    case 'vk'%In Vetterli and Kovacevic book
        if di==1;h0=[1 2 1]/4;h1=[-1 -2 6 -2 -1]/4;
        else h0=[-1 2 6 2 -1]/4;h1=[-1 2 -1]/4;end	        
        h0=ftrans2(h0,t);h1=ftrans2(h1,t);% McClellan transfrom to obtain 2D diamond filters
    case {'ko','kos'}%Orthogonal filters in Kovacevic's thesis (non smooth and smooth versions)
        if strcmp(DFtype,'ko');a0=2;a1=0.5;a2=1;else a0=-sqrt(3);a1=-sqrt(3);a2=2+sqrt(3);end        	
        h0=[ 0  -a1       -a0*a1 0;
            -a2 -a0*a2    -a0    1;
             0   a0*a1*a2 -a1*a2 0];	
        % h1 = qmf2(h0);
        h1=[0 -a1*a2 -a0*a1*a2 0;
            1  a0    -a0*a2    a2;
            0 -a0*a1  a1       0];	       
        norm=sqrt(2)/sum(h0(:)); % Normalize filter sum and norm;	
        h0=h0*norm;h1=h1*norm;	      
        if ~di%Reverse filters for reconstruction
            for n=1:2;h0=flip(h0,n);h1=flip(h1,n);end
        end    
    case 'lax'%By Lu, Antoniou and Xu
        h=[-1.2972901e-5  1.2316237e-4 -7.5212207e-5  6.3686104e-5  9.4800610e-5 -7.5862919e-5  2.9586164e-4 -1.8430337e-4;
            1.2355540e-4 -1.2780882e-4 -1.9663685e-5 -4.5956538e-5 -6.5195193e-4 -2.4722942e-4 -2.1538331e-5 -7.0882131e-4;	    
           -7.5319075e-5 -1.9350810e-5 -7.1947086e-4  1.2295412e-3  5.7411214e-4  4.4705422e-4  1.9623554e-3  3.3596717e-4;	    
	        6.3400249e-5 -2.4947178e-4  4.4905711e-4 -4.1053629e-3 -2.8588307e-3  4.3782726e-3 -3.1690509e-3 -3.4371484e-3;	    
            9.6404973e-5 -4.6116254e-5  1.2371871e-3 -1.1675575e-2  1.6173911e-2 -4.1197559e-3  4.4911165e-3  1.1635130e-2;	    
           -7.6955555e-5 -6.5618379e-4  5.7752252e-4  1.6211426e-2  2.1310378e-2 -2.8712621e-3 -4.8422645e-2 -5.9246338e-3;	    
            2.9802986e-4 -2.1365364e-5  1.9701350e-3  4.5047673e-3 -4.8489158e-2 -3.1809526e-3 -2.9406153e-2  1.8993868e-1;	    
           -1.8556637e-4 -7.1279432e-4  3.3839195e-4  1.1662001e-2 -5.9398223e-3 -3.4467920e-3  1.9006499e-1  5.7235228e-1];	
        h0=sqrt(2)*[h               h(:,end-1:-1:1);
                    h(end-1:-1:1,:) h(end-1:-1:1,end-1:-1:1)];		
        h1=mirrorModulation(h0);	
    case 'sk'%By Shah and Kalker
        h=[ 0.621729    0.161889   -0.0126949  -0.00542504  0.00124838;
	        0.161889   -0.0353769  -0.0162751  -0.00499353  0;
	       -0.0126949  -0.0162751   0.00749029  0           0;
           -0.00542504  0.00499353  0           0           0;
            0.00124838  0           0           0           0];	
        h0=sqrt(2)*[h(end:-1:2,end:-1:2) h(end:-1:2,:);
                    h(:,end:-1:2)        h];	
        h1=mirrorModulation(h0);	
    case 'dvmlp'
        q=sqrt(2);b=.02;b1=b*b;
        h=[b/q  0        -2*q*b 0        3*q*b 0       -2*q*b  0        b/q;
           0   -1/(16*q)  0     9/(16*q) 1/q   9/(16*q)  0    -1/(16*q) 0;
           b/q  0        -2*q*b 0        3*q*b 0       -2*q*b  0        b/q];     
        g0=[-b1/q 0        4*b1*q            0          -14*q*b1            0           28*q*b1             0          -35*q*b1             0           28*q*b1             0          -14*q*b1            0           4*b1*q           0       -b1/q;
             0    b/(8*q)  0                -13*b/(8*q)  b/q                33*b/(8*q) -2*q*b              -21*b/(8*q)  3*q*b              -21*b/(8*q)  -2*q*b              33*b/(8*q)  b/q               -13*b/(8*q)  0                b/(8*q)  0;
            -q*b1 0       -1/(256*q)+8*q*b1  0           9/(128*q)-28*q*b1 -1/(q*16)   -63/(256*q)+56*q*b1  9/(16*q)    87/(64*q)-70*q*b1   9/(16*q)   -63/(256*q)+56*q*b1 -1/(q*16)    9/(128*q)-28*q*b1  0          -1/(256*q)+8*q*b1 0       -q*b1;
             0    b/(8*q)  0                -13*b/(8*q)  b/q                33*b/(8*q) -2*q*b              -21*b/(8*q)  3*q*b              -21*b/(8*q) -2*q*b               33*b/(8*q)  b/q               -13*b/(8*q)  0                b/(8*q)  0;
            -b1/q 0        4*b1*q            0          -14*q*b1            0           28*q*b1             0          -35*q*b1             0           28*q*b1             0          -14*q*b1            0           4*b1*q           0       -b1/q];
        if di;h1=mirrorModulation(g0);h0=h;
        else h1=mirrorModulation(h);h0=g0; 
        end
    case {'cd', '7-9'}%By Cohen and Daubechies
        % 1D prototype filters: the '7-9' pair
        h0=[0.026748757411 -0.016864118443 -0.078223266529 0.266864118443 0.602949018236 0.266864118443 -0.078223266529 -0.016864118443 0.026748757411];
        g0=[-0.045635881557 -0.028771763114 0.295635881557 0.557543526229 0.295635881557 -0.028771763114 -0.045635881557];
        if di;h1=mirrorModulation(g0,2);
        else h1=mirrorModulation(h0,2);h0=g0;
        end        
        h0=sqrt(2)*ftrans2(h0,t);h1=sqrt(2)*ftrans2(h1,t);%Use McClellan to obtain 2D filters
    case 'sinc'%The "sinc" case, no Perfect Reconstruction. Ideal low and high pass filters
        flength = 30;	
        h0=fir1(flength,0.5);
        h1=mirrorModulation(h0,2);	
        h0=sqrt(2)*ftrans2(h0,t);h1=sqrt(2)*ftrans2(h1,t);%Use McClellan to obtain 2D filters        
        if ~di% Reverse filters for reconstruction
            for n=1:2;h0=flip(h0,n);h1=flip(h1,n);end
        end     
    case 'oqf_362'%Some "home-made" filters!
        h0=sqrt(2)/64*[ sqrt(15)   -3   0;
                        0           5  -sqrt(15);
                       -2*sqrt(15)  30  0;
                        0           30  2*sqrt(15);
                        sqrt(15)    5   0;
                        0          -3   -sqrt(15)]';
         h1=-mirrorModulation(h0);
         for n=1:2;h1=flip(h1,n);end%h1 = -reverse2(modulate2(h0, 'b'));, I've interpreted as before, it was commented
         if ~di%Reverse filters for reconstruction
             for n=1:2;h0=flip(h0,n);h1=flip(h1,n);end
         end
    case 'test'%Only for the shape, not for PR
        h0=[0 1 0;1 4 1;0 1 0];h1=[0 -1 0;-1 4 -1;0 -1 0];	
    case 'testDVM'%Only for directional vanishing moment
        h0=[1 1;1 1]/sqrt(2);h1=[-1 1;1 -1]/sqrt(2);		 	
    case 'qmf'%By Lu, Antoniou and Xu (ideal response / window)
        %Window
        m=2;n=2;
        w=[];h=[];
        w1d=kaiser(4*m+1,2.6);
        for n1=-m:m
            for n2=-n:n
                w(n1+m+1,n2+n+1)=w1d(2*m+n1+n2+1)*w1d(2*m+n1-n2+1);
                h(n1+m+1,n2+n+1)=.5*sinc((n1+n2)/2)*.5*sinc((n1-n2)/2);
            end
        end      
        h=sqrt(2)*h/sum(h(:));
        h0=h.*w;h1=mirrorModulation(h0);    
     case 'qmf2'%By Lu, Antoniou and Xu (ideal response / window)    
        h=[-.001104  .002494 -.001744  .004895 -.000048 -.000311;
            .008918 -.002844 -.025197 -.017135  .003905 -.000081;
           -.007587 -.065904  .100431 -.055878  .007023  .001504;
            .001725  .184162  .632115  .099414 -.027006 -.001110;
           -.017935 -.000491  .191397 -.001787 -.010587  .002060;
            .001353  .005635 -.001231 -.009052 -.002668  .000596];
        h0=h./sum(sum(h));h1=mirrorModulation(h0);
    case {'dmaxflat4','dmaxflat5','dmaxflat6','dmaxflat7'}
        M1=1/sqrt(2);M2=M1;
        k1=1-sqrt(2);k3=k1;k2=M1;          
        h=[.25*k2*k3 .5*k2 1+.5*k2*k3]*M1;h=[h fliplr(h(1:end-1))];
        g=[-.125*k1*k2*k3 0.25*k1*k2 -0.5*k1-0.5*k3-0.375*k1*k2*k3 1+.5*k1*k2]*M2;g=[g fliplr(g(1:end-1))];
        B=dmaxflat(str2double(DFtype(end)),0);
        h0=mctrans(h,B);g0=mctrans(g,B);
        h0=sqrt(2)*(h0./sum(h0(:)));g0=sqrt(2)*(g0./sum(g0(:)));       
        h1=mirrorModulation(g0);
        if di==0;h1=mirrorModulation(h0);h0=g0;end   
    otherwise%Assume the "degenerated" case: 1D wavelet filters	
        if di==1;[h0,h1]=wfilters(DFtype,'d');else [h0,h1]=wfilters(DFtype,'r');end
end

if modu;h0=mirrorModulation(h0/sqrt(2),2);h1=[];end

function h = dmaxflat(N,d)
% returns 2-D diamond maxflat filters of order 'N' 
% the filters are nonseparable and 'd' is the (0,0) coefficient, being 1 or 0 depending on use
% by Arthur L. da Cunha, University of Illinois Urbana-Champaign
% Aug 2004
% Part of the Nonsubsampled Contourlet Toolbox
% (http://www.mathworks.de/matlabcentral/fileexchange/10049-nonsubsampled-contourlet-toolbox)

assert(N>=1 && N<=7,'N must be in {1,2,3,4,5,6,7}');
if N==1;h=[0 1 0;1 0 1;0 1 0]/4;
elseif N==2;h=[0 -1 0;-1 0 10;0 10 0]/32;
elseif N==3;h=[0 3 0 2;3 0 -27 0;0 -27 0 174;2 0 174 0]/512;
elseif N==4;h=[0 -5 0 -3 0; -5 0 52 0 34;0 52 0 -276 0;-3 0 -276 0 1454;0 34 0 1454 0]/2^12;
elseif N==5;h=[0 35 0 20 0 18;35 0 -425 0 -250 0;0 -425 0 2500 0 1610;20 0 2500 0 -10200 0;0 -250 0 -10200 0 47780;18 0 1610 0 47780 0]/2^17;
elseif N==6;h=[0 -63 0 -35 0 -30 0;-63 0 882 0 495 0 444;0 882 0 -5910 0 -3420 0;-35 0 -5910 0 25875 0 16460;0 495 0 25875 0 -89730 0;-30 0 -3420 0 -89730 0 389112;0 44 0 16460 0 389112 0]/2^20;
elseif N==7;h =[0 231 0 126 0 105 0 100;231 0 -3675 0 -2009 0 -1715 0;0 -3675 0 27930 0 15435 0 13804;126 0 27930 0 -136514 0 -77910 0;0 -2009 0 -136514 0 495145 0 311780;105 0 15435 0 495145 0 -1535709 0;0 -1715 0 -77910 0 -1535709 0 6305740;100 0 13804 0 311780 0 6305740 0]/2^24;
end
if N>1;h=[h fliplr(h(:,1:end-1))];h=[h;flipud(h(1:end-1,:))];end
h(N+1,N+1)=d;

end

function h = mctrans(b,t)
% MCTRANS McClellan transformation
%   H = mctrans(B,T) produces the 2-D FIR filter H that
%   corresponds to the 1-D FIR filter B using the transform T.
% Convert the 1-D filter b to SUM_n a(n) cos(wn) form
% Part of the Nonsubsampled Contourlet Toolbox
% (http://www.mathworks.de/matlabcentral/fileexchange/10049-nonsubsampled-contourlet-toolbox)

nn=(length(b)-1)/2;
b=rot90(fftshift(rot90(b,2)),2); % Inverse fftshift
a=[b(1) 2*b(2:nn+1)];
inset=floor((size(t)-1)/2);

% Use Chebyshev polynomials to compute h
P0=1;P1=t;
h=a(2)*P1;
rows=inset(1)+1;cols=inset(2)+1;
h(rows,cols)=h(rows,cols)+a(1)*P0;
for l=3:nn+1
    P2=2*conv2(t,P1);
    rows=rows+inset(1);cols=cols+inset(2);
    P2(rows,cols)=P2(rows,cols)-P0;
    rows=inset(1)+(1:size(P1,1));cols=inset(2)+(1:size(P1,2));
    hh=h;
    h=a(l)*P2;h(rows,cols)=h(rows,cols)+hh;
    P0=P1;P1=P2;
end
h = rot90(h,2);%Rotate for use with filter2

end

end