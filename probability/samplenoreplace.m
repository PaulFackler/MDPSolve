% samplenoreplace Draws samples without replacement
% USAGE
%   I=samplenoreplace(m,n,p);
% INPUTS
%   m  : # of elements in set
%   n  : # of samples
%   p  : # of elements in ach sample (p<=m)
% OUTPUT
%   I  : p x n mattrix of values on {1,...,m}
% Each column contains p unique values 
function I=samplenoreplace(m,n,p)
e=rand(p,n);

if exist('samplenoreplacec','file')
  I=samplenoreplacec(e,m);
  return
end

I=zeros(p,n);
ind=(1:m);
for j=1:n
  for i=0:p-1
    i1=i+1;
    k=ceil((m-i)*e(i1,j));
    ik=ind(k);
    ind(k)=ind(m-i);
    ind(m-i)=ik;
    I(i1,j)=ik;
  end
end
return