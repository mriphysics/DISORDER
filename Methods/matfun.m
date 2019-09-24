function A=matfun(f,A,B)

%MATFUN   Overrides the behaviour of pagefun for non-gpu arrays
%   A=MATFUN(F,A,{B})
%   * F is a function handle. Interesting handles are @ctranspose, @inv,
%   @mldivide, @mrdivide, @mtimes, @transpose
%   * A is the first input array
%   * {B} is the second input array
%   ** A is the output array
%

gpu=isa(A,'gpuArray');

if (isequal(f,@mldivide) || isequal(f,@mrdivide)) && gpu && verLessThan('matlab','R2016a')%%%EXTEND THIS TO OTHER METHODS
    f=@inv;
    if isequal(f,@mrdivide);[B,A]=parUnaFun({A,B});end
end

NDA=numDims(A);
if NDA<=2
    if nargin<3
        A=f(A);return
    else
        NDB=numDims(B);
        if NDB<=2;A=f(A,B);return;end
    end
end
perm=1:NDA;perm([2 1])=1:2;

if isa(A,'gpuArray')      
    if nargin>=3        
        if isequal(f,@mtimes)
            NA=size(A);NB=size(B);
            if NA(1)==1 && NB(2)==1
                A=permute(A,perm);                           
                A=sum(bsxfun(@times,A,B),1);
                if NB(2)==1;A=permute(A,perm);end
            else
                A=pagefun(f,A,B);
            end
        else
            A=pagefun(f,A,B);
        end
    else
        A=pagefun(f,A);
    end
else
    if isequal(f,@ctranspose) || isequal(f,@transpose)        
        A=permute(A,perm);
        if isequal(f,@ctranspose);A=conj(A);end
    elseif ~isempty(A) 
        isB=(nargin>=3 && ~isempty(B));
        if isB
            [NA,NB]=parUnaFun({A,B},@size);
            [NDA,NDB]=parUnaFun({A,B},@numDims);
            [NDA,NDB]=parUnaFun({NDA,NDB},@max,3);
            NA(end+1:NDB)=1;
            NB(end+1:NDA)=1;
            NAB=[1 1 NA(3:end)./NB(3:end)];
            NBA=[1 1 NB(3:end)./NA(3:end)];
            B=repmat(B,ceil(NAB));
            A=repmat(A,ceil(NBA));
        end
        NA=size(A);
        [A,NAP]=resSub(A,3:numDims(A));NAP(end+1:3)=1;        
        if ~isB            
            for p=1:NAP(3);A(:,:,p)=f(A(:,:,p));end
        else
            B=resSub(B,3:numDims(B));
            C=zeros([NA(1) NB(2) NAP(3)],'like',A);
            for p=1:NAP(3);C(:,:,p)=f(A(:,:,p),B(:,:,p));end
            A=C;
        end
        A=resSub(A,3,NA(3:end));
    end
end
