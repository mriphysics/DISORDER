function x=fold(x,m,NS,NY,A,fl)

%FOLD   SENSE folding operation along a given dimension
%   X=FOLD(X,M,NS,NY,{A},{FL})
%   * X is the array on which to operate
%   * M is the direction across which to operate
%   * NS is the spatial size of the array
%   * NY is the spectral size of the array
%   * {A} is a matrix to perform folding
%   * {FL} is a flag to implement mirrored folding. Not tested that it 
%   works well when the "undersampling" is bigger than 3, but it should be 
%   of little interest for the main application at the moment which is to
%   implement mirror boundary conditions
%   ** X is the folded array
%

if nargin<5;A=[];end
if nargin<6 || isempty(fl);fl=0;end

if NS~=NY   
    ND=ndims(x);
    oddFactSENSE=2*ceil((ceil(NS/NY)-1)/2)+1;
    oFRed=(oddFactSENSE-3)/2;    
    oFRedT=oFRed*NY;
    over=(NS-NY)/2;
    FOV=[floor(over) ceil(over)];
    if mod(NS,NY)~=0 || fl%Slower
        if ~isempty(A) && ~fl
            x=aplGPU(A,x,m);
        else   
            if ND<=14                
                if m==1;xb=x(FOV(2)+1:NS-FOV(1),:,:,:,:,:,:,:,:,:,:,:,:,:);
                elseif m==2;xb=x(:,FOV(2)+1:NS-FOV(1),:,:,:,:,:,:,:,:,:,:,:,:);
                elseif m==3;xb=x(:,:,FOV(2)+1:NS-FOV(1),:,:,:,:,:,:,:,:,:,:,:);
                else xb=dynInd(x,FOV(2)+1:NS-FOV(1),m);
                end
            else
                xb=dynInd(x,FOV(2)+1:NS-FOV(1),m);
            end
            for s=1:oFRed
                if mod(s,2)==1 && fl
                    if ND<=14
                        if m==1;xb=xb+flip(x(FOV(2)+1-s*NY:NS-FOV(1)-s*NY,:,:,:,:,:,:,:,:,:,:,:,:,:)+x(FOV(2)+1+s*NY:NS-FOV(1)+s*NY,:,:,:,:,:,:,:,:,:,:,:,:,:),m);
                        elseif m==2;xb=xb+flip(x(:,FOV(2)+1-s*NY:NS-FOV(1)-s*NY,:,:,:,:,:,:,:,:,:,:,:,:)+x(:,FOV(2)+1+s*NY:NS-FOV(1)+s*NY,:,:,:,:,:,:,:,:,:,:,:,:),m);
                        elseif m==3;xb=xb+flip(x(:,:,FOV(2)+1-s*NY:NS-FOV(1)-s*NY,:,:,:,:,:,:,:,:,:,:,:)+x(:,:,FOV(2)+1+s*NY:NS-FOV(1)+s*NY,:,:,:,:,:,:,:,:,:,:,:),m);
                        else xb=xb+flip(dynInd(x,FOV(2)+1-s*NY:NS-FOV(1)-s*NY,m)+dynInd(x,FOV(2)+1+s*NY:NS-FOV(1)+s*NY,m),m);
                        end
                    else
                        xb=xb+flip(dynInd(x,FOV(2)+1-s*NY:NS-FOV(1)-s*NY,m)+dynInd(x,FOV(2)+1+s*NY:NS-FOV(1)+s*NY,m),m);
                    end
                else
                    if ND<=14
                        if m==1;xb=xb+x(FOV(2)+1-s*NY:NS-FOV(1)-s*NY,:,:,:,:,:,:,:,:,:,:,:,:,:)+x(FOV(2)+1+s*NY:NS-FOV(1)+s*NY,:,:,:,:,:,:,:,:,:,:,:,:,:);
                        elseif m==2;xb=xb+x(:,FOV(2)+1-s*NY:NS-FOV(1)-s*NY,:,:,:,:,:,:,:,:,:,:,:,:)+x(:,FOV(2)+1+s*NY:NS-FOV(1)+s*NY,:,:,:,:,:,:,:,:,:,:,:,:);
                        elseif m==3;xb=xb+x(:,:,FOV(2)+1-s*NY:NS-FOV(1)-s*NY,:,:,:,:,:,:,:,:,:,:,:)+x(:,:,FOV(2)+1+s*NY:NS-FOV(1)+s*NY,:,:,:,:,:,:,:,:,:,:,:);
                        else xb=xb+dynInd(x,FOV(2)+1-s*NY:NS-FOV(1)-s*NY,m)+dynInd(x,FOV(2)+1+s*NY:NS-FOV(1)+s*NY,m);
                        end
                    else
                        xb=xb+dynInd(x,FOV(2)+1-s*NY:NS-FOV(1)-s*NY,m)+dynInd(x,FOV(2)+1+s*NY:NS-FOV(1)+s*NY,m);
                    end
                end
            end            
            if mod(oFRed,2)==0 && fl%Flip may only work for even differences between NS and NY 
                if ND<=14
                    if m==1
                        xb(1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:,:)=xb(1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:,:)+flip(x(1:FOV(2)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:,:),m);
                        xb(NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:,:,:)=xb(NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:,:,:)+flip(x(NS-FOV(1)+1+oFRedT:NS,:,:,:,:,:,:,:,:,:,:,:,:,:),m);
                    elseif m==2
                        xb(:,1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:)=xb(:,1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:)+flip(x(:,1:FOV(2)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:),m);
                        xb(:,NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:,:)=xb(:,NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:,:)+flip(x(:,NS-FOV(1)+1+oFRedT:NS,:,:,:,:,:,:,:,:,:,:,:,:),m);
                    elseif m==3
                        xb(:,:,1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:)=xb(:,:,1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:)+flip(x(:,:,1:FOV(2)-oFRedT,:,:,:,:,:,:,:,:,:,:,:),m);
                        xb(:,:,NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:)=xb(:,:,NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:)+flip(x(:,:,NS-FOV(1)+1+oFRedT:NS,:,:,:,:,:,:,:,:,:,:,:),m);
                    else
                        xb=dynInd(xb,1:FOV(1)-oFRedT,m,dynInd(xb,1:FOV(1)-oFRedT,m)+flip(dynInd(x,1:FOV(2)-oFRedT,m),m));
                        xb=dynInd(xb,NY-FOV(2)+1+oFRedT:NY,m,dynInd(xb,NY-FOV(2)+1+oFRedT:NY,m)+flip(dynInd(x,NS-FOV(1)+1+oFRedT:NS,m),m));
                    end                       
                else                    
                    xb=dynInd(xb,1:FOV(1)-oFRedT,m,dynInd(xb,1:FOV(1)-oFRedT,m)+flip(dynInd(x,1:FOV(2)-oFRedT,m),m));
                    xb=dynInd(xb,NY-FOV(2)+1+oFRedT:NY,m,dynInd(xb,NY-FOV(2)+1+oFRedT:NY,m)+flip(dynInd(x,NS-FOV(1)+1+oFRedT:NS,m),m));
                end
            else
                if ND<=14
                    if m==1
                        xb(1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:,:)=xb(1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:,:)+x(NS-FOV(1)+1+oFRedT:NS,:,:,:,:,:,:,:,:,:,:,:,:,:);
                        xb(NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:,:,:)=xb(NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:,:,:)+x(1:FOV(2)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
                    elseif m==2
                        xb(:,1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:)=xb(:,1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:)+x(:,NS-FOV(1)+1+oFRedT:NS,:,:,:,:,:,:,:,:,:,:,:,:);
                        xb(:,NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:,:)=xb(:,NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:,:)+x(:,1:FOV(2)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:,:);                        
                    elseif m==3
                        xb(:,:,1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:)=xb(:,:,1:FOV(1)-oFRedT,:,:,:,:,:,:,:,:,:,:,:)+x(:,:,NS-FOV(1)+1+oFRedT:NS,:,:,:,:,:,:,:,:,:,:,:);
                        xb(:,:,NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:)=xb(:,:,NY-FOV(2)+1+oFRedT:NY,:,:,:,:,:,:,:,:,:,:,:)+x(:,:,1:FOV(2)-oFRedT,:,:,:,:,:,:,:,:,:,:,:,:);                        
                    else
                        xb=dynInd(xb,1:FOV(1)-oFRedT,m,dynInd(xb,1:FOV(1)-oFRedT,m)+dynInd(x,NS-FOV(1)+1+oFRedT:NS,m));
                        xb=dynInd(xb,NY-FOV(2)+1+oFRedT:NY,m,dynInd(xb,NY-FOV(2)+1+oFRedT:NY,m)+dynInd(x,1:FOV(2)-oFRedT,m));
                    end
                else
                    xb=dynInd(xb,1:FOV(1)-oFRedT,m,dynInd(xb,1:FOV(1)-oFRedT,m)+dynInd(x,NS-FOV(1)+1+oFRedT:NS,m));
                    xb=dynInd(xb,NY-FOV(2)+1+oFRedT:NY,m,dynInd(xb,NY-FOV(2)+1+oFRedT:NY,m)+dynInd(x,1:FOV(2)-oFRedT,m));
                end
            end
            x=xb;
        end
    else
        N=size(x);        
        x=reshape(x,[N(1:m-1) NY ceil(NS/NY) prod(N(m+1:ND))]);
        x=sum(x,m+1);
        x=reshape(x,[N(1:m-1) NY N(m+1:ND)]);
        perm=zeros(1,ND);perm(m)=-FOV(2);     
        x=circshift(x,perm);
    end
end
