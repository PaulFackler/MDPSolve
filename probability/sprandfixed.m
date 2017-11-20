% sprandfixed Sparse random matrix with fixed number of non-zeros per column or row
% USAGE
%   A=sprandfixed(m,n,p,rows);
% INPUTS
%   m,n : size of output
%   p   : # of non-zeros per column
%   rows : 0/1 if 1 fixed # of non-zeros per row [default: 0]
% OUTPUTS
%   A   : m x n matrix
% 
% non-zero elements of A are distributed uniformly on [0,1]
function A=sprandfixed(m,n,p,rows)
if nargin<4 || isempty(rows), rows=false; end
if rows
  if p>n, error('p must be <= n'); end
  I=samplenoreplace(n,m,p);
  A=sparse(repmat(1:m,p,1),I,rand(m,p),m,n);
else
  if p>m, error('p must be <= m'); end
  I=samplenoreplace(m,n,p);
  A=sparse(I,repmat(1:n,p,1),rand(p,n),m,n);
end


  
