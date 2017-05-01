% getnetworkP cpomputes a transition probability matrix for a network of variables
% The children variables must be conditionallly independent given the parents
% USAGE
%   [P,X,Ix,S] = getnetworkP(n,children,parents,p);
% INPUTS
%  n        : N vector with sizes of each variables 
%               (may include variables not in current network)
%  children : nc vector of variable indices on {1,...,N}
%  parents  : nc cell array of row vectors variable indices on {1,...,N}
%  p        : nc element cell array with kth element 
%                n(children(k)) x prod(n(parents{children(k)}))
%  order    : nc vector, a permutation of integers 1 to nc indicating the order
%               that children variables should be processed
% OUTPUT
%   P       : R x C transition matrix for Prob(children|parents) - R=prod(n(children))
%   X       : C x np matrix of possible combinations of parent values  
%   Ix      : C-elment index vector for children values in X
%   S       : R x nc matrix of possible combinations of children values  
%
% It is assumed that the rows of P correspond to the children in sorted order
% i.e., sort(children) and the columns correspond to the parents in sorted order, 
% i.e., pv=unique([parents{:}]) where np=length(pv) and C=prod(n(pv))
%
% See also: getnetworkEV
function [P,X,Ix,S]=getnetworkP(n,children,parents,p)
if exist('order','var') && ~isempty(order) 
  children=children(order); parents=parents(order); p=p(order);
end
nc=length(children);
include=unique([parents{:}]);
ni=length(include);
if nargout<3
  X=networkgetX(n,children,parents);
else
  [X,Ix,S]=networkgetX(n,children,parents);
end
ii=zeros(1,length(n)); ii(include)=1:ni;
P=p{1}(:,getI(X,ii(parents{1})));
for i=2:nc
  P=kroncol(P,p{i}(:,getI(X,ii(parents{i}))));
end
return


