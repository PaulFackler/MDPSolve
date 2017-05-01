% Dynamic programming algorithm to determine optimal chain operation
% With d objects chain operations require d-1 binary operations
% and the cost of an operation, given by cost(i,j,k), is the cost of operating
% on an object that has merged objects i through j-1 with an object that has merged 
% objects j through k (note: i<j<=k)
% This is appropriate when multiple intermerdiate products, e.g. ((A,B),(C,D))
% are allowed.
%
% Given a d-element cell array A with A{i} a matrix of size m(i) x n(i)
% i.e., [m,n]=cellfun(@size,A);
% Use with
% for i=1:d-1, A{order(i,1)}=operator(A{order(i,1)},A{order(i,2)}); end
% A{1} will contain the result (note this will overwrite the A{i}).
% For chain matrix multiplication: cost=@(i,j,k) m(i)*m(j)*n(k)
% For Kronecker products: cost=@(i,j,k) prod(m(i:k))*prod(n(i:k))
% For columnwise Kronecker (Khatri-Rao) products: cost=@(i,j,k) prod(m(i:k))*n
%
% These cost functions assume that the A{i} are full. If they are sparse the cost
% of Kronecker products is prod(q(i:k)) where q=cellfun(@nnz,A); 
% If we assume that there is at least 1 non-zero per row and column then for 
% chain matrix multiplication the following cost function could be used
%    function c=mccost(i,j,k,m,n,q)
%       if i==j-1
%          if j==k
%            c=q(i)*q(j);
%          else
%            c=q(i)*n(k);
%          end
%       else
%          if j==k
%            c=m(i)*q(j);
%          else
%            c=m(i)*m(j)*n(k);
%          end
%       end
% With Khatri-Rao the cost depends on the pattern of sparcity but bounds or 
% expectations based on random placement could be computed.
function [order,c]=optchainord(d,cost)
% initialize C by putting pairwise costs on the super-diagonal
C=zeros(d,d);
for i=1:d
  for j=i+1:d
    C(i,j)=cost(i,j,j);
  end
end
S=diag(2:d,1);               % best split for (i,j) (split is first element in second group) 
for len=3:d
  i1=len;
  for i0=1:d-len+1
    [c,s]=min(C(i0,i0:i1-1)+C(i0+1:i1,i1)');
    C(i0,i1)=c+cost(i0,s,i1);
    S(i0,i1)=i0+s;
    i1=i1+1;
  end
end

% get sequence
order=zeros(0,2);
order=getorder(1,d,order,S);
c=C(1,d); % total cost of sequence
return

% recursive function to recover optimal order
function order=getorder(F,L,order,S)
  SFL=S(F,L);
  if SFL>0, 
    if F<SFL-1, order=getorder(F,SFL-1,order,S);end
    if SFL<L,   order=getorder(SFL,L,order,S);  end
    order=[order;F SFL];
  end