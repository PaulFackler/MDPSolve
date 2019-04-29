% rowselect Selects the value of a specified row for each column of a matrix
% USAGE
%   [x,ind]=rowselect(X,ind);
% INPUTS
%   X   : m x n matrix
%   ind : n-vector of integers in {1,...,m}
% OUTPUT
%   x   : 1 x n vector with x(j) = X(ind(j),j)
%   ind : 1 x n vector of linear index values
function [x,ind]=rowselect(X,ind)
[rX,cX] = size(X);
ind = ind(:)' + (0:rX:rX*(cX-1));
x=X(ind);