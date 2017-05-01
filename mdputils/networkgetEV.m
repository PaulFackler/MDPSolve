% getnetworkEv sets up an expected value function for a network of variables
% The children variables must be conditionallly independent given the parents
% USAGE
%   EV=getnetworkEV(n,children,parents,p,order);
% INPUTS
%   n        : N vector with sizes of each variables 
%                (may include variables not in current network)
%   children : nc vector of variable indices on {1,...,N}
%   parents  : nc cell array of row vectors variable indices on {1,...,N}
%   p        : nc element cell array with kth element 
%                n(children(k)) x prod(n(parents{children(k)}))
%                Currently all p{i} must be full (not sparse)
%   order    : nc vector, a permutation of integers 1 to nc indicating the order
%                that children variables should be processed
% OUTPUT
%   EV      : EV function where EV(v) returns an m-vector equal to E[v(children)|parents]
%               m equals prod(n(I)) where I is the union of all the parent vectors
%
% It is assumed that v is in sorted order, i.e., sort(children)
% and the output is in sorted order, i.e., unique([parents{:}])
%
% !!!! NOTE: this may fail if any of the n are equal to 1 !!!!!
%
% USES tprod: https://www.mathworks.com/matlabcentral/fileexchange/16275-tprod-arbitary-tensor-products-between-n-d-arrays
function EV=getnetworkEV(n,children,parents,p,order)
if exist('order','var') && ~isempty(order) 
  children=children(order); parents=parents(order); p=p(order);
end
n=n(:)'; children=children(:)'; % row vectors
N=length(n);
nc=length(children);
rind=sort(children);
EV=getevfunc(N,nc,n,rind,children,parents,p);
return

% get function handle for EV function
function EV=getevfunc(N,nc,n,rind,children,parents,p)
EV=@(v)computeev(v,N,nc,n,rind,children,parents,p);

% this is the actual function called when EV(v) is invoked
function ev=computeev(v,N,nc,n,rind,children,parents,p)
ev=v;
cvars=false(1,N);    % logical vector cvars(i)=1 if variable i is currently in output 
for k=1:nc
  child=children(k);
  ev=calltrod(ev,p{k},[n n],rind,find(cvars),child,parents{k},child);
  rind(rind==child)=[];   % input variables included  
  cvars(parents{k})=true;
  ev=squeeze(ev);
end
ev=ev(:);
return

function AB=calltrod(A,B,n,rvarsa,cvarsa,rvarsb,cvarsb,sumvars)
N=length(n);
aind=[fliplr(rvarsa) fliplr(cvarsa)];
bind=[fliplr(rvarsb) fliplr(cvarsb)];
A=reshape(A,n(aind));
B=reshape(B,n(bind));
N1=N+1;
N21=2*N+1;
aind=[N1-fliplr(rvarsa) N21-fliplr(cvarsa)];
bind=[N1-fliplr(rvarsb) N21-fliplr(cvarsb)];

[~,ii]=ismember(N1-sumvars,aind);
if length(ii)~=length(sumvars)
  error('sumvars not correctly specified')
end
aind(ii)=-aind(ii);
[~,ii]=ismember(N1-sumvars,bind);
if length(ii)~=length(sumvars)
  error('sumvars not corrently specified')
end
bind(ii)=-bind(ii);
AB=tprod(A,aind,B,bind);
