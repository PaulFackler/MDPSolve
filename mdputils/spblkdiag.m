% spblkdiag Sparse block diagonal matrix
% USAGE 
%   X=spblkdiag(x1,x2,...);
% or
%   X=spblkdiag(C);
% INPUTS
%   x1, x2, ..., xn  : n numeric matrices
%   C            : cell array of n numeric matrices = {x1,x2,...,xn}
% OUTPUT
%   X            : block diagonal matrix 
%
%  X = [x1  0 ...  0;
%        0 x2 ...  0;
%             ...
%        0  0 ... xn];
                                                  
function X=spblkdiag(x1,varargin)
if nargin<1
  X=[];
else
  if isnumeric(x1)
    x1=sparse(x1);
    X=[{x1} varargin];
  elseif iscell(x1)
    x1{1}=sparse(x1{1});
    X=[x1 varargin];
  else
    error('inproper inputs');
  end
    X=blkdiag(X{:});
end
  