function ind=electrostaticRepulsionDISORDER(N)

%ELECTROSTATICREPULSIONDISORDER computes a permutation of a given 2D grid
%using electrostatic repulsion for sequencing
%   IND=ELECTROSTATICREPULSIONDISORDER(N)
%   * N is the size of the grid
%   ** IND is the computed ordering
%

poi=zeros(prod(N),2);
ind=zeros(prod(N),1);

indFree=ones(1,prod(N));

ind(1)=sub2ind(N,ceil(N(1)/2),ceil(N(2)/2));
[poi(1,1),poi(1,2)]=ind2sub(N,ind(1));

indFree(ind(1))=0;
Nshift{1}(:,1)=[-N(1);0;N(1)];
Nshift{2}(1,:)=[-N(2) 0 N(2)];

for n=2:prod(N)
    poiD=reshape(poi(1:n-1,:),[1 1 n-1 2]);
    poiD=repmat(poiD,[3 3 1 1]);
    for l=1:2
        poiD(:,:,:,l)=bsxfun(@plus,poiD(:,:,:,l),Nshift{l});
    end
    poiD=reshape(poiD,[9*(n-1) 2]);
    %Minimum of the wrapped electric repulsion
    minRepul=inf;
    indMinRepul=[];
    for m=1:prod(N)
        if indFree(m)
            [poiCand(1,1),poiCand(1,2)]=ind2sub(N,m);
            dist=sum(bsxfun(@minus,poiCand,poiD).^2,2);      
            wdist=min(reshape(dist,[9 n-1]),[],1);
            wrepul=sum(1./wdist);
            if wrepul<minRepul-(1e-9)%To avoid numerical issues
                minRepul=wrepul;
                indMinRepul=m;
            end
        end
    end
    ind(n)=indMinRepul;
    [poi(n,1),poi(n,2)]=ind2sub(N,ind(n));
    indFree(ind(n))=0;
end





            