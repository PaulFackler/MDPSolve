% evalnetwork evaluates the expected value function for a network of variables
% The children variables must be conditionallly independent given the parents
% USAGE
%   ev=evalnetwork(v,n,children,parents,p);
% INPUTS
%  v        : prod(n(children)) vector
%  n        : N vector with sizes of each variables 
%               (may include variables not in current network)
%  children : nc vector of variable indices on {1,...,N}
%  parents  : nc cell array of vectors variable indices on {1,...,N}
%  p        : nc element cell array with kth element 
%                n(children(k)) x prod(n(parents{children(k)}))
% OUTPUT
%   ev      : m-vector equal to E[v(children)]parents]
%               m equals prod(n(I)) where I is the union of all the parent vectors
%   sizes   : size of the output array at the end of each iteration 
%               (can be used to evaluate the amount of work performed)
%
% This can be used to create an EV function for MDPSolve using
%   EV = @(v) evalnetwork(v,p,children,parents,n);
%
% It is assumed that v is in sorted order, i.e., sort(children)
% and the output is in sorted order, i.e., unique([parents{:}])
function [ev,sizes]=evalnetwork(v,n,children,parents,p)
n=n(:)'; children=children(:)'; % row vectors
N=length(n);
maxel=numel(v);
ev=v;
nc=length(children);
rind=sort(children);
cvars=false(1,N);    % logical vector cvars(i)=1 if variable i is currently in output 
sizes=zeros(nc,1);
for k=1:nc
  child=children(k);
  ev=calltrod(ev,p{k},[n n],rind,find(cvars),child,parents{k},child);
  rind(rind==child)=[];   % input variables included  
  cvars(parents{k})=true;
  ev=squeeze(ev);
  if maxel<numel(ev), maxel=numel(ev); end
  sizes(k)=numel(ev);
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


function AB=callbsxfun(A,B,n,rvarsa,cvarsa,rvarsb,cvarsb,sumvars)
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