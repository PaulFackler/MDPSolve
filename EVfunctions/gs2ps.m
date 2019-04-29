% gs2ps Converts a factored model defined by functions to one defined by CPTs
% USAGE
%   [p,pparents]=gs2ps(g,gparents,s,X,e,options)
% INPUTS
%   g        : m-element cell array of transition functions
%   gparents : m-element cell array of conditioning variable (parent) indices
%   s        : cell array of state variables
%   X        : matrix or cell array of state/action variables
%   e        : cell array of rv structures (discrete or w/ discrete approximations)
%   options  : structure variable (fields described below)
% OUTPUTS
%   p        : m-element cell array of conditional probability matrices
%   pparents : m-element cell array of conditioning variable (parent) indices
%
% parent vectors use positive numbers to refer to the X variables and
%   negative numbers to refer to the e variables
% gparents and pparents may differ because the variable order for the CPTs
%   is always ordered first by the X variables and then by the e variables
%   and because any e variables associated with only a single state
%   variable are summed out
%
% Options fields:
%   cleanup  : method used to handle extrapolation (see g2P)
%
% 
function [p,pparents]=gs2ps(g,gparents,s,X,e,options)
% initializations
if nargin<6, options = struct(); end
if nargin<5, e={}; end
if isa(g,'function_handle'); g={g}; end
if isnumeric(gparents); gparents={gparents}; end
if isnumeric(s); s={s}; end
if isstruct(e); e={e}; end

cleanup = 0;
if isfield(options,'cleanup'), cleanup = options.cleanup; end

% get variable dimensions
ds = length(g);
if isnumeric(X);
  dx=size(X,2);
else
  dx = length(X);
end
de = length(e);

% determine which noise variables can be summed out
A=zeros(ds,de);
for i=1:ds
  pxi=gparents{i}(gparents{i}>0);
  if any(pxi<1 | pxi>dx)
    error(['parents{' num2str(i) '} contains incorrect values']);
  end
  pei=-gparents{i}(gparents{i}<0);
  if any(pei<1 | pei>de)
    error(['parents{' num2str(i) '} contains incorrect values']);
  end
  A(i,pei)=1;
  if ~isa(g{i},'function_handle')
    error('g elements must be function handles')
  end
  g{i} = rvdef('f',g{i},s{i});
end
summable = sum(A)<=1;
p=cell(1,ds);
pparents=cell(1,ds);

useD=false;
if useD
  % convert factored model to diagram
  D=f2d(g,gparents,s,X,e);

  % call conditional to compute the CPTs
  for i=1:ds
    pxi=gparents{i}(gparents{i}>0);
    % get list of noise variables that cannot be summed
    pei = false(1,de);
    pei(-gparents{i}(gparents{i}<0)) = true;
    pei(summable) = false;
    pei = find(pei);
    p{i}=conditional(D,dx+de+i,[pxi dx+pei],options);
    pparents{i} = [pxi -pei];
  end
else % convert to CPTs directly 
  evals    = cellfun(@(e) e.values, e, 'UniformOutput',false);
  eweights = cellfun(@(e) e.cpt,    e, 'UniformOutput',false);
  ne = cellfun(@length,evals);
  for i=1:ds
    pxi=gparents{i}(gparents{i}>0);
    % get list of noise variables that cannot be summed
    pei = false(1,de);
    pei(-gparents{i}(gparents{i}<0)) = true;
    xevals = cell(1,length(pxi) + sum(pei));
    if iscell(X)
      [xevals{:}]=rectgrid(X{pxi}, evals{pei});
    else
      [xevals{:}]=rectgrid(unique(X(:,pxi),'rows'), evals(pei));
    end
    gvals=double(g{i}.valfunc(xevals{:}));
    p{i} = rectbas(gvals,s{i},[],cleanup);
    sume = find(pei(summable));
    ww=eweights(sume);
    sume2=matchindices(sume,find(pei))+1;
    p{i} = sumout(p{i},[size(p{i},2)/prod(ne(pei)) ne(pei)],sume2,ww);
    pei(summable) = false;
    pei = find(pei);
    pparents{i} = [pxi -pei];
  end
end
