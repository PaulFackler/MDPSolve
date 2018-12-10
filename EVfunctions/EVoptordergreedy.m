% greedy method of determining optimal order and grouping
% compares gamma =(N-n)/(N*n) where N is the size of the new variables
% added and n is the size of the future state eliminated. The variable with 
% the smallest gamma value is added

% tis should be adjusted to account for e variables that get summed out.
function [order]=EVoptordergreedy(p,parents,X,e,options)
factorsize = true;
if nargin>=5 && isstruct(options)
  if isfield(options,'factorsize'), factorsize = options.factorsize; end
end
m=length(parents);
ns=cellfun(@(x)size(x,1),p);
if iscell(X)
  nx=cellfun(@(x)length(x),X);
end
ne=cellfun(@(x)length(x.values),e);

A = zeros(m,length(e));
for i=1:m
  A(-parents{i}(parents{i}<0),i)=1;
end
numused = sum(A); % number of times e(j) appears

notprocessed = 1:m;  % indices of the remaining variables
numleft = m;         % number of variables not yet processed (length(notprocessed))
common = [];         % list of conditioning variables already entered
order = zeros(1,m);
for i=1:m-1
  mingam=inf;
  for j=1:numleft
    ii = notprocessed(j);
    % remove already entered variables
    jj = matchindices(common,parents{ii});
    if ~isempty(jj)
      parents{ii}(jj) = []; 
    end
    % size of new X variables 
    if iscell(X)
      N = prod(nx(parents{ii}(parents{ii}>0))) ;
    else
      [~,Xi] = getI(X,parents{ii}(parents{ii}>0));
      N=size(Xi,1);
    end
    % size of new e variables; if factorsize include only ones not summed
    % out at this step
    evars = -parents{ii}(parents{ii}<0);
    if factorsize && any(numused(evars)==1)
      evars(matchindices(find(numused(evars)==1),evars)) = [];
    end
    N = N * prod(ne(evars));
    % gamma value of current variable
    gam =  (N-ns(ii))/ (N*ns(ii));
    if gam < mingam
      mingam = gam;
      minv = ii;
    end
  end
  order(i) = minv;
  common = union(common,parents{minv});
  notprocessed(matchindices(minv,notprocessed)) = [];
  numleft = numleft - 1;
  evars = -parents{minv}(parents{minv}<0);
  numused(evars) = numused(evars) - 1;
end
order(m) = notprocessed;
      
 