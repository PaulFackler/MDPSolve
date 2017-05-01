% bottom up O(d^2) algorithm
% note this is only optimal if we allow only one intermediate product
% because we examine whether it is better to add to the front end or the back end
% but does not allow, e.g. ((A,B),(C,D))
function [ind,c]=chainkroncol(m)
d=length(m);
% initialize C by putting pairwise costs on the super-diagonal
C=diag(m(1:d-1).*m(2:d),1); % C(i,j) is the minimum cost of linking i-j
E=false(d,d);               % true if end is added, false if start is added 
for len=3:d
  i1=len;
  for i0=1:d-len+1
    m01=prod(m(i0:i1));
    c0=m01+C(i0+1,i1); % cost of adding first last
    c1=m01+C(i0,i1-1); % cost of adding last last
    if c0>c1 % add last last
      E(i0,i1)=true; 
      C(i0,i1)=c1;
    else     % add first last
      C(i0,i1)=c0;
    end
    i1=i1+1;
  end
end

% get sequence
ind=zeros(1,d);
i=1; j=d;
for k=d:-1:3
  if E(i,j)
    ind(k)=j;
    j=j-1;
  else
    ind(k)=i;
    i=i+1;
  end    
end
ind(1)=i;
ind(2)=j;
c=C(1,d); % total cost of sequence
return