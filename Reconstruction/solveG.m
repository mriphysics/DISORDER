function g=solveG(E,EH,A,C,R,nIt,tolType,tol,deb,estGF,y,x)

%SOLVEG   Computes the g-factors for a given reconstruction
%   G=SOLVEG(E,EH,{A},{C},{R},{NIT},{TOLTYPE},{TOL},{DEB},{ESTGF},{Y},{X})
%   * E is an encoding structure
%   * EH is a decoding structure
%   * {A} is a preconditioner structure
%   * {C} is a constrain structure
%   * {R} is a regularization structure
%   * {NIT} is the maximum number of iterations
%   * {TOLTYPE} is the type of tolerance used for convergence
%   * {TOL} is the tolerance set for convergence
%   * {DEB} indicates whether to print information about convergence
%   * {ESTGF} indicates the mode of estimation
%   * {Y} is the measured data
%   * {X} is the reconstructed data
%   ** G is a g-factor map
%

if nargin<3;A=[];end
if nargin<4;C=[];end
if nargin<5;R=[];end
if nargin<6 || isempty(nIt);nIt=300;end
if nargin<7 || isempty(tolType);tolType='Energy';end%tolType='NormwiseBackward2Error'
if nargin<8 || isempty(tol);tol=1e-2;end
if nargin<9 || isempty(deb);deb=2;end
if nargin<10;estGF=[1 0];end
if nargin<11;x=[];end
if nargin<12;y=[];end

rmc=0;
if ~isfield(C,'Ma')
    C.Ma=single(sum(abs(E.Sf),4)~=0);
    C.Ma=dynInd(C.Ma,2:size(C.Ma,6),6,0);
    rmc=1;
end
g=C.Ma;
Caux=C;
if rmc;Caux=[];end
NG=size(g);

if estGF(1)<=1%Analytic computation
    indPEFull=cell(1,3);
    for n=1:3
        indPE=[];
        oddFactSENSE=2*ceil((ceil(E.Uf{n}.NX/E.Uf{n}.NY)-1)/2)+1;
        oFRed=(oddFactSENSE-3)/2;    
        oFRedT=oFRed*E.Uf{n}.NY;
        over=(E.Uf{n}.NX-E.Uf{n}.NY)/2;
        FOV=[floor(over) ceil(over)];
        if n==3 && isfield(E,'Bf')%MultiBand
            MB=E.Uf{n}.NX/E.Uf{n}.NY;
            for m=1:MB;indPE{m}=(m-1)*E.Uf{n}.NY+1:m*E.Uf{n}.NY;end
        elseif n==E.pe && (isfield(E,'Bf') && estGF(1)==0.5) %&& (mod(E.Uf{n}.NY,(2*E.Bl))~=0 || any(C.Ma(:)~=1))%Necessary to split for MultiBand, otherwise g-factor maps with ringing- If estGF~=0.5 we assume there's not much ringing or it can be prevented by other means...
           for m=1:E.Uf{n}.NX;indPE{m}=m;end
        else        
            if n~=E.pe || ~(isfield(E,'Bf') && estGF(1)==0.5)
                if oddFactSENSE>1
                    cont=1;
                    indPE{cont}=1:FOV(2)-oFRedT;
                    cont=cont+1;
                    for s=-oFRed:oFRed
                        indPE{cont}=FOV(2)+1+E.Uf{n}.NY*s:FOV(2)+E.Uf{n}.NY*(s+1);
                        cont=cont+1;
                    end
                    indPE{cont}=E.Uf{n}.NX-FOV(1)+1+oFRedT:E.Uf{n}.NX;
                    cont=cont+1;
                else
                    indPE{1}=1:NG(n);
                end
            else
                if oddFactSENSE>1
                    cont=1;
                    aux=1:FOV(2)-oFRedT;             
                    NF=E.Uf{n}.NY/(2*E.Bl);         
                    assignPE;
                    for s=-oFRed:oFRed
                        aux=FOV(2)+1+E.Uf{n}.NY*s:FOV(2)+E.Uf{n}.NY*(s+1);                 
                        assignPE;
                    end
                    aux=E.Uf{n}.NX-FOV(1)+1+oFRedT:E.Uf{n}.NX;                
                    assignPE;             
                else
                    aux=1:E.Uf{n}.NX;
                    assignPE;
                end
            end                
        end
        indPEFull{n}=indPE;
    end
end

yG=g;
if estGF(1)>1%Monte-Carlo computation
    yR=yG;yR(:)=0;
    for r=1:ceil(estGF(1)-1)
        xGF=abs(CGsolver(y+plugNoise(y)/sqrt(2),E,EH,A,Caux,R,[],nIt,tolType,tol,0,0,deb)-x).^2;
        yR=yR+xGF;
    end    
    yR=dynInd(yR,1,6);
    wF=ones(1,3)*estGF(2);wF(end+1:numDims(x))=1;wF=wF(1:numDims(x));
    wC=ones(wF,'like',yR);
    yR=padarray(yR,double((wF-1)/2),'replicate');
    yR=convn(yR,wC,'valid');
    g=dynInd(g,1,6).*sqrt((yR/(ceil(estGF(1)-1)*prod(wF))).*sum(abs(dynInd(E.Sf,1,6)).^2,4));
else
    for n=1:length(indPEFull{1})
        for m=1:length(indPEFull{2})        
            for o=1:length(indPEFull{3})
                indIt={indPEFull{1}{n},indPEFull{2}{m},indPEFull{3}{o},1};
                Wp=dynInd(C.Ma,indIt,[1:3 6]);          
                if any(Wp(:)~=0)             
                    yG(:)=0;yG=dynInd(yG,indIt,[1:3 6],1);                  
                    gaux=encode(yG,E);
                    gaux=constrain(decode(gaux,EH,E),C);
                    gaux=CGsolver(gaux,E,EH,A,Caux,R,[],nIt,tolType,tol,-1,0,deb);              
                    g=dynInd(g,indIt,[1:3 6],dynInd(gaux,indIt,[1:3 6]));                
                end
            end
        end
    end    
    g=dynInd(g,1,6);
    g=sqrt(abs(g));
end

function assignPE
    ind=1;
    while ind<=length(aux)
        indPE{cont}=aux(ind:min(ind+NF-1,length(aux)));
        ind=ind+NF;
        cont=cont+1;
    end
end

end

