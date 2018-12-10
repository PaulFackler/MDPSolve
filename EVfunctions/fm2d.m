% f2d Converts a factored into an influence diagram specification
% USAGE
%    D=fm2d(p,parents,s,X,e,fnames,xnames,enames)
% INPUTS
%   p        : m-element cell array of conditional probability matrices
%   parents  : m-element cell array of conditioning variable (parent) indices
%   X        : matrix or cell array of state/action variables
%   e        : cell array of rv structures (discrete or w/ discrete approximations)
%   fnames   : cell array of names for the target ('f') variables [optional]
%   xnames   : cell array of names for the conditioning ('X') variables [optional]
%   enames   : cell array of names for the noise ('e') variables [optional]
%                If the names are not passed variables will be called
%                f1, f2, ..., x1, x2, ... , e1, e2 ,... 
% OUTPUT
%   D        : a diagram structure (see add2diagram for description)
%
% The diagram structure created can be passed as is to either
%   conditional or condexp to obtain a transition matrix or EV function.
% The following modifications are required before it is passed to d2model:
%   all conditioning variables are state ('s') type; action variables need 
%     to be altered to 'a' type
%   future state ('f') type variables should be linked to current state ('s')
%     variables by giving them the same names with a '+' appended
%   a reward/utility variable must be added to D
%
% The EV functions created by EVcreate and condexp use different computational
%   approachs but will given results that differ only by rounding error  
% Similarly the transition matrices created by fm2P and conditional use
%   different computational approachs but will given results that differ 
%   only by rounding error
function D=fm2d(p,parents,s,X,e,fnames,xnames,enames)
if ~iscell(e), e={e}; end
% get variable dimensions
ds=length(p);
if isempty(s), s=cell(1,ds); end
dx=size(X,2);
de=length(e);
% set default names
if nargin<6
  fnames=cell(1,ds);
  for i=1:ds
    fnames{i} = ['f' num2str(i)];
  end
  xnames=cell(1,dx);
  for i=1:dx
    xnames{i} = ['x' num2str(i)];
  end
  enames=cell(1,de);
  for i=1:de
    enames{i} = ['e' num2str(i)];
  end
end

if iscell(X)
  nx = cellfun(@length,X);
else
  nx = zeros(1,dx);
  for i=1:dx
    nx(i) = length(unique(X(:,i)));
  end
end

% create diagram stucture ordering variable types as x, e, f
D=[];
for i=1:dx
  if iscell(X)
    D=add2diagram(D,xnames{i},'s',1,{},X{i});
  elseif ~isempty(X)
    D=add2diagram(D,xnames{i},'s',1,{},unique(X(:,i)));
  else
    D=add2diagram(D,xnames{i},'s',1,{},(1:nx(i))');
  end
end

for i=1:de
  D=add2diagram(D,enames{i},'c',0,{},e{i});
end

for i=1:ds
  mi = length(parents{i});
  q=cell(1,mi);
  for j=1:mi
    if parents{i}(j) > 0, qj = parents{i}(j);
    else                  qj = dx-parents{i}(j);
    end
    q{j}=D.names{qj};
  end
  if isnumeric(p{i})
    if isempty(s{i})
      s{i}=(1:size(p{i},1))';
    end
    cpt = rvdef('d',p{i},s{i});
  elseif isa(p{i},'function_handle')
    cpt = rvdef('f',p{i},s{i});
  else
    cpt = p{i};
  end
  D=add2diagram(D,fnames{i},'f',1,q,cpt);
end