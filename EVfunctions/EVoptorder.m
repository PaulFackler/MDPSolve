% EVoptorder Finds the optimal order to combine factors & processes a network
% USAGE
%   [order,mergevec]=EVoptorder(P,parents,X,n)
% INPUTS
%   P        : ds cell array of CPTs
%   parents  : ds element cell array of index vectors of conditioning variables
%   X        : nx x dx matrix of conditioning variable values
%   n        : dx vector of sizes of conditioning variables (used to determine optimal order)
%                If omitted conditioning variables sizes will be determined from X
% OUTPUTS
%   order    : variable order 
%   mergevec : grouping vector 
function [order,mergevec]=EVoptorder(p,parents,X,e,options)
torder=[];
usegreedy=true;
if exist('options','var')
  if isfield(options,'torder'), torder=options.torder; end
  if isfield(options,'usegreedy'), usegreedy=options.usegreedy; end
end
if isnumeric(X)
  error('X must be a cell array to determine the optimal order')
end
ds=length(p);         % # of state variables
dx=size(X,2);         % # of conditioning variables
de=length(e);         % # of noise variables

% get variable sizes
n=[cellfun(@(x)length(x),x) cellfun(@(x)size(x,1),p)];

% get the optimal processing order
if isempty(torder)
  V=zeros(ds+1,dx+ds);  % vertex matrix
  V(1,dx+1:dx+ds)=1;
  for i=1:ds
    V(i+1,parents{i})=1; 
    V(i+1,dx+i)=1;
  end
  Q=[ones(1,dx) zeros(1,ds)];
  options=struct('printlevel',0,'alg',1);
  torder=ndsumprodorder(V,Q,n,options)-1;
end

ps=num2cell(1:ds);  % cell array to hold reordering information
order=[];
for i=1:ds
  o1=torder(i,1);
  o2=torder(i,2);
  if o1==0, 
    order=[order ps{o2}]; %#ok<AGROW>
  else
    ps{o1}=[ps{o1} ps{o2}];
    ps{o2}=[];
  end
end
mergevec=cellfun(@length,ps(order));
mergevec=mergevec(mergevec>0);