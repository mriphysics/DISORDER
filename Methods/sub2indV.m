function ind=sub2indV(N,sub)

%SUB2INDV Implements a vectorial version of the sub2ind function
%   IND=SUB2INDV(N,SUB)
%   * N is the size of the multidimensional array the subscripts correspond to
%   * SUB is a matrix of subscripts arranged as # indexes times # dimensions
%   ** IND is the resulting set of indexes
%

sub=double(sub);
assert(size(sub,2)==numel(N),'Number of dimensions of subscripts (%d) does not match dimensions sizes (%d)',size(sub,2),numel(N));
N=N(:)';
sublim=bsxfun(@minus,sub,N);
assert(gather(all(sub(:)>0 & sublim(:)<=0)),'Subscripts outside the range given by dimensions sizes');
Ncum=[1 cumprod(N(1:end-1))];
ind=sum(bsxfun(@times,(sub-1),Ncum),2)+1;
