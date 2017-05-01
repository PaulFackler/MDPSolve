% sprandC sparse random matrix with set # of non-zeros per column or row
% USAGE
%   A=sprandC(m,n,p,rows);
% INPUTS
%   m    : # of rows
%   n    : # of columns
%   p    : # of non-zeros per column
%   rows : 1 for p non-zeros per row [optional]
% OUTPUT
%   A : m x n sparse matrix with nnz(A)=n*p (or m*p is rows=1)

function A=sprandC(m,n,p,rows)
if nargin<4, rows=false; end
if rows
  ii=repmat(1:m,p,1);
  jj=zeros(m,p);
  for i=1:m
    jj(i,:)=randperm(n,p);
  end
  A=sparse(ii,jj,rand(m,p),m,n);
else
  jj=repmat(1:n,p,1);
  ii=zeros(p,n);
  for i=1:n
    ii(:,i)=randperm(m,p);
  end
  A=sparse(ii,jj,rand(p,n),m,n);
end
