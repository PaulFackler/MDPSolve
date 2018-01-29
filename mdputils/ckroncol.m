% ckroncol Repeated columnwise Kronecker products on a cell array of matrices.
% Solves (B1xB2x...xBd) where x denotes columnwise Kronecker (Khatri-Rao) product.
% USAGE
%   B=ckroncol(A)   
% INPUT 
%   A   :  a d-element cell array of matrices with each element a matrix
%            with n columns
% OUTPUT
%   B   : a prod(m) x n matrix where m is a d-vector containing the row
%            numbers of the A matrices

function z=ckroncol(B)
[m,n]=cellfun(@size,B);
if any(n~=n(1))
  error('B matrices must all have the same # of columns')
end
z=B{1};
for i=2:length(m)
  z=kroncol(z,B{i});
end