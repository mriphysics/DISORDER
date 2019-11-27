function x=ifold(x,m,NS,NY,A,fl)

%IFOLD   Inverse SENSE folding operation along a given dimension
%   X=IFOLD(X,M,NS,NY,{A})
%   * X is the array on which to operate
%   * M is the direction across which to operate
%   * NS is the spatial size of the array
%   * NY is the spectral size of the array
%   * {A} is a matrix to perform inverse folding
%   * {FL} is a flag to implement inverse mirrored folding
%   ** X is the unfolded array
%

if nargin<5;A=[];end
if nargin<6 || isempty(fl);fl=0;end

if NS~=NY
    ND=ndims(x);
    oddFactSENSE=2*ceil((ceil(NS/NY)-1)/2)+1;
    oFRed=oddFactSENSE-2;
    disc=(oddFactSENSE*NY-NS)/2;
    iFOV=[floor(disc) ceil(disc)]; 
    
    if mod(NS,NY)~=0 || fl~=0%Slower
        if ~isempty(A) && fl==0
            x=aplGPU(A,x,m);
        else
            if fl==0
                indUnf=[1+iFOV(1):NY repmat(1:NY,[1 oFRed]) 1:NY-iFOV(2)];
            elseif mod((oFRed-1)/2,2)==0
                indUnf=repmat(1:NY,[oFRed 1]);
                indUnf(2:2:end,:)=flip(indUnf(2:2:end,:),2);
                indUnf=indUnf';     
                indUnf=[flip(1:NY-iFOV(2)) indUnf(:)' flip(1+iFOV(1):NY)];
            else
                indUnf=repmat(1:NY,[oFRed 1]);
                indUnf(1:2:end,:)=flip(indUnf(1:2:end,:),2);
                indUnf=indUnf';     
                indUnf=[1+iFOV(1):NY indUnf(:)' 1:NY-iFOV(2)];
            end
            if ND<=14
                if m==1;x=x(indUnf,:,:,:,:,:,:,:,:,:,:,:,:,:);
                elseif m==2;x=x(:,indUnf,:,:,:,:,:,:,:,:,:,:,:,:);
                elseif m==3;x=x(:,:,indUnf,:,:,:,:,:,:,:,:,:,:,:);
                else x=dynInd(x,indUnf,m);
                end
            else
                x=dynInd(x,indUnf,m);
            end
        end
    else
        perm=zeros(1,ND);perm(m)=-iFOV(1);
        x=circshift(x,perm);
        perm=ones(1,ND);perm(m)=NS/NY;
        x=repmat(x,perm);                
    end
end
