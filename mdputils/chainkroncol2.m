% bottom up O(d^3) algorithm
% allows multiple intermerdiate products, e.g. ((A,B),(C,D))
% use with
% for i=1:d-1, A{order(i,1)}=kroncol(A{order(i,1)},A{order(i,2)}); end
% A{1} will contain the result
function [order,c]=chainkroncol(m)
d=length(m);
% initialize C by putting pairwise costs on the super-diagonal
C=diag(m(1:d-1).*m(2:d),1)+diag(m,0); % C(i,j) is the minimum cost of linking i-j
S=diag(2:d,1);               % best split for (i,j) (split is first element in second group) 
for len=3:d
  i1=len;
  for i0=1:d-len+1
    [c,s]=min(C(i0,i0:i1-1)+C(i0+1:i1,i1)');
    C(i0,i1)=c+prod(m(i0:i1));
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