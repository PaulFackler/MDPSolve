% amdpweights computes AM weights from historical data
% USAGE
%   w=amdpweights(P,X,Xhist,svars,w0);
% INPUTS
%   P     : q-element cell array of alternative transition matrices
%   X     : matrix of discrete X values or cell array of vectors if X is rectanglar
%   Xhist : T-row matrix of historical values of X 
%   svars : vector  of columns of X for state variables
%   w0    : initial weighting vector (1 x q)
% OUTPUT
%   w     : Txq matrix of historical weights (note: w0 is the first row)
function w=amdpweights(P,X,Xhist,svars,w0)
T=size(Xhist,1);
q=numel(w0);
w=zeros(T,q);
% rectangular grid
if iscell(X)
  Bs=rectbas(Xhist(2:T,svars),X(svars),[],2);
  Bx=rectbas(Xhist(1:T-1,:),X,[],2);
% X is composed of arbitrary points
else
  Bs=scatterbas(Xhist(2:T,svars),tesselate(X(:,svars)));
  Bx=scatterbas(Xhist(1:T-1,:),tesselate(X));
end
w(1,:)=w0(:)';
W=zeros(T-1,q);
for j=1:q
   W(:,j)=sum((P{j}*Bx).*Bs,1)';
end
for t=1:T-1
  wt=W(t,:).*w(t,:);
  wt=wt/sum(wt);
  w(t+1,:)=wt;
end