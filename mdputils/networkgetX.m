% getnetworkX computes the matrix of all parent values for a network of variables
% USAGE
%   X = getnetworkX(n,children,parents);
% INPUTS
%   n        : N vector with sizes of each variables 
%                (may include variables not in current network)
%   children : nc vector of variable indices on {1,...,N}
%   parents  : nc cell array of row vectors variable indices on {1,...,N}
% OUTPUT
%   X       : rX x cX matrix of possible combinations of parent values  
%
% The parents are listed in the sorted order defined by n, i.e., pv=unique([parents{:}]) 
%   rX=prod(n(pv)) 
%   cX=length(pv) 
%
% See also: getnetworkEV, getnetworkP
function [X,Ix,S]=getnetworkX(n,children,parents)
include=unique([parents{:}]);
cX=length(include);
X=cell(1,cX);
for i=1:cX
  X{i}=(1:n(include(i)))';
end
X=rectgrid(X);

if nargout>1
  [~,cols]=ismember(children,include);
  [Ix,S]=getI(X,cols);
end