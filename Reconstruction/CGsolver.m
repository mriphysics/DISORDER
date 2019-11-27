function [x,En,gF,cF,nn,rF]=CGsolver(y,E,EH,A,C,R,x,nIt,tolType,tol,estGF,estCF,deb)

%CGSOLVER   Performs a CG-based pseudoinverse reconstruction
%   [X,EN,GF,CF,NN]=CGSOLVER(Y,E,EH,{A},{C},{R},{X},{GF},{NIT},{TOLTYPE},{TOL},{ESTGF},{ESTCF},{DEB},{SINGLECH})
%   * Y is the array of measures
%   * E is an encoding structure
%   * EH is a decoding structure
%   * {A} is a preconditioner structure
%   * {C} is a constrain structure
%   * {R} is a regularizer structure
%   * {X} is an initial condition for the array to be reconstructed
%   * {NIT} is the maximum number of iterations
%   * {TOLTYPE} is the type of tolerance used for convergence
%   * {TOL} is the tolerance set for convergence
%   * {ESTGF} indicates whether to estimate the g-factor
%   * {ESTCF} indicates whether to estimate the chi2-factor
%   * {DEB} indicates whether to print information about convergence
%   ** X is the reconstruction
%   ** EN are the energies of successive reconstructions
%   ** GF is a g-factor value
%   ** CF is the chi2 factor
%   ** NN is the number of iterations for convergence
%   ** RF are the residuals
%

if nargin<4;A=[];end
if nargin<5;C=[];end
if nargin<6;R=[];end
if nargin<8 || isempty(nIt);nIt=300;end%300;end
if nargin<9 || isempty(tolType);tolType='Energy';end%tolType='NormwiseBackward2Error'
if nargin<10 || isempty(tol);tol=1e-5;end%5e-3;end%1e-2;end
if nargin<11 || isempty(estGF);estGF=[0 0];end
if nargin<12 || isempty(estCF);estCF=0;end
if nargin<13 || isempty(deb);deb=1;end
cF=[];gF=[];rF=[];

NitTest=2;
tolTyp=stopCondition(tolType,'CG');

if isfield(A,'Ps')%Optimal circular preconditioner from RH Chan, "Circulant preconditioners for Hermitian Toeplitz systems," SIAM J Matrix Anal Appl, 10(4), 542-550, 1989
    if isfield(R,'Fd');A.Ps=1./(1+mean(R.Fd.la(:))*conj(A.Ps).*A.Ps);end
    if isfield(R,'Fo');A.Ps=1./(1+mean(R.Fo.la(:))*conj(A.Ps).*A.Ps);end        
end

if estGF(1)~=-1
    b=constrain(decode(y,EH,E),C);
    if isfield(R,'Tp') && isfield(R.Tp,'x0');b=b+regularize(x,R,3);end
else
    b=constrain(y,C);
    if tolTyp==6;tolTyp=0;end
end

if nargin<7 || isempty(x);x=b;x(:)=0;EHE=x;
elseif ~isfield(E,'Je') || E.Je==0;EHE=constrain(decode(encode(x,E),EH,E)+regularize(x,R),C);
else EHE=constrain(encodedecode(x,E,EH)+regularize(x,R),C);
end

if ~ismember(tolTyp,6:7)
    zb=precondition(b,A);%For normalizing the stopping condition
    zb=conj(zb).*b;%For normalizing the stopping condition
    bnorm2=sum(abs(zb(:)));%For normalizing the stopping condition
end

r=b-EHE;b=[];%We can enter y-Ex to be decoded and something empty in x to save computation time 
z=precondition(r,A);
p=z;
zr=conj(z).*r;
rsold=sum(zr(:));
En=[];
if all(abs(rsold(:))<1e-6)
    if estCF;cF=x;end
    nn=0;
    fprintf('Early termination of CG, data probably corrupted!\n');
    return;    
end
err=inf;
if tolTyp==6;En=inf(2,nIt);end

for n=1:max(nIt)    
    if ~isfield(E,'Je') || E.Je==0;EHE=constrain(decode(encode(p,E),EH,E)+regularize(p,R),C);else EHE=constrain(encodedecode(p,E,EH)+regularize(p,R),C);end
    
    g=rsold./sum(conj(p(:)).*EHE(:));%Previous implementation was using conj(rsold) here 
    xup=bsxfun(@times,g,p);
    x=x+xup;
    
    if mod(n,NitTest)==1 && tolTyp==6
        En(:,n)=computeEnergy(y,x,E,R,EH);
        if deb>=2;fprintf('It %d - Ene Fid: %0.2g / Ene Reg: %0.2g / Ene Tot: %0.2g\n',n,En(1,n),En(2,n),sum(En(:,n)));end
    end
    r=r-bsxfun(@times,g,EHE);
    z=precondition(r,A);
    zr=conj(z).*r;

    rs=sum(zr(:));
    d=rs./(rsold+eps);
    p=z+bsxfun(@times,d,p);
    rsold=rs;

    if tolTyp==1;err=sqrt(sum(abs(zr(:)))/bnorm2);
    elseif tolTyp==6 && mod(n,NitTest)==1 && n>1;err=sum(En(:,n-NitTest)-En(:,n))/sum(En(:,n-NitTest));%err=sqrt(max(0,err));
    end
    if deb>=2 && (tolTyp<6 || mod(n,NitTest)==1);fprintf('It %d - Err %s: %0.2g\n',n,tolType,err);end
    if err<tol || all(abs(rs(:))<1e-6);break;end 
end
nn=n;

if n==max(nIt) && deb==2;fprintf('CG solver terminated without reaching convergence\n');end
if tolTyp==6;En(:,(1:nIt)>n)=[];end
En=En(:,1:NitTest:end);

%COMPUTE THE CHI2-FACTORS
if estCF;[cF,rF]=solveC(x,y,E,EH);end

%COMPUTE THE G-FACTORS
if estGF(1)>0
    if deb>=1;fprintf('Number of iterations for g-factor calculations: %d\n',nn);end
    gF=solveG(E,EH,A,C,R,nn,tolType,0,deb,estGF,y,x);
end

