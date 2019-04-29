% EVcreate Creates an EV function from conditional probabilities
% USAGE
%   [EV,workspace]=EVcreate(p,parents,X,e,options);
% INPUTS
%   p       : m-element cell array of conditional probability matrices
%   parents : m-element cell array of conditioning variable (parent) indices
%   X       : matrix or cell array of state/action variables
%   e       : cell array of rv structures (discrete or w/ discrete approximations)
%   options : structure variable (fields described below)
% OUTPUTS
%   EV        : a function handle for an EV function
%   workspace : a structure variable containing the workspace of the EV function
%
% parents{i} contains the columns of X that condition the ith variable
% Thus Xi=unique(X(:,parents{i}),'rows') contains the unique elements of the
%   conditioning variables that condition output variable i.
% It must be the case that the number of rows in Xi matches the number 
%   of columns in p{i}: size(Xi,1)=size(p{i},2)
%
% Options fields:
%   order       : m-vector if alternative variable order is desired
%   getoptord   : attempts to find optimal combination method
%   mergevec    : vector providing merge information with sum(mergevec)=m
%                   An EV function is created with length(mergevec) components
%                   The default is ones(1,m), i.e., one component for each variable
%                   mergevec=m results in a single transition matrix
%                   Pass as empty to avoid using the optimal grouping (same
%                     as a vector of 1s)
%
%   useI        : if 0 uses J indexing from the beginning
%   alg         : the algorithm used by indexedmult 
%   vreorder    : reordering vector for v (ignored if order option is defined)
%   nx          : sizes of the X variables (used to determine optimal grouping)
%   expandall   : expands all of the CPTS to have nx columns
%
% Note: in general the workspace output is not needed but is provided here for  
% debugging purposes or to help curious users understand how this function works
%
% If EVcreate is called multiple times using different values of p but with
%   no change to X, parents or options, then call it the first time using
%     [EV,workspace]=EVcreate(p,parents,X,e,options);
% and on subsequent calls use
%     EV=EVcreate(p,workspace);
function [EV,ws]=EVcreate(p,parents,X,e,options)
if nargin<4, e={}; end
if nargin<5, options=[]; end
order=[];
getoptord=false;
getoptgroup=true;
mergevec=ones(1,length(p));  % default is no merging
pm=[];
nx=[];
expandall=false;
alg=[];
useI=true;
vreorder=[];
ws=[];
Ionly = false; 
if exist('options','var') && ~isempty(options)
  if isfield(options,'getoptord'),   getoptord   = options.getoptord;   end
  if isfield(options,'mergevec'),    mergevec    = options.mergevec;    end
  if isfield(options,'pm'),          pm          = options.pm;          end
  if isfield(options,'nx'),          nx          = options.nx;          end
  if isfield(options,'expandall'),   expandall   = options.expandall;   end
  if isfield(options,'vreorder'),    vreorder    = options.vreorder;    end
  if isfield(options,'order'),       order       = options.order;       end
  if isfield(options,'usebsxfun'),   alg         = options.alg;         end
  if isfield(options,'useI'),        useI        = options.useI;        end
  if isfield(options,'Ionly'),       Ionly       = options.Ionly;       end
  if isfield(options,'workspace'),   ws          = options.workspace;   end
end

% if mergevec defined in options as empty then then no grouping will
% be performed
if exist('mergevec','var')
  if isempty(mergevec), getoptgroup=false; end
end

m=length(parents);
% check if mergevec is correctly specified
if ~isempty(mergevec) && sum(mergevec)~=m,
  error('mergevec is incorrectly specified - must sum to m')
end

alg=3;
% use this to force factor 2 to use J
useI = true;
firstJ=1;
% J not working for models with noise
%useI=true;
%firstJ = ifthenelse(isempty(e),0,m+1);  
ws = fprocess(p,parents,X,e,options);
p=ws.p;
w=ws.w;
Ip=ws.Ip;
Iy=ws.Iy;
Jp=ws.Jp;
Jy=ws.Jy;
vreorder = ws.vreorder;
EV=@(varargin) EVeval(p,w,Ip,Iy,Jp,Jy,vreorder,useI,alg,firstJ,varargin{:}); 

% packworkspace Creates a workspace structure
function workspace=packworkspace(p,parents,e,w,Ip,Iy,Jp,Jy,vreorder,mergevec)
  workspace.p = p;
  workspace.parents = parents;
  workspace.e = e;
  workspace.w = w;
  workspace.Ip=Ip;
  workspace.Iy=Iy;
  workspace.Jp=Jp;
  workspace.Jy=Jy;
  workspace.vreorder=vreorder;
  workspace.mergevec=mergevec;
  
% unpackworkspace Extracts the elements of a workspace structure
function [w,Ip,Iy,Jp,Jy,vreorder,mergevec]=unpackworkspace(workspace)
   w=workspace.w;
  Ip=workspace.Ip;
  Iy=workspace.Iy;
  Jp=workspace.Jp;
  Jy=workspace.Jy;
  vreorder=workspace.vreorder;
  mergevec=workspace.mergevec;
  
% fprocess Gets the variables needed to form an EVfunction from a factored model
% USAGE
%   ws =fprocess((p,parents,X,e);
% INPUTS
%   p       : m-element cell array of conditional probability matrices
%   parents : m-element cell array of conditioning variable (parent) indices
%   X       : matrix or cell array of state/action variables
%   e       : cell array of rv structures (discrete or w/ discrete approximations)
% OUTPUTS
%   ws      : workspace structure containing the following fields:
%               p    - m-element cell array of CPTs
%               Ip   - m-element cell array of expansion indices for p
%               Iy   - m-element cell array of expansion indices for p
%               Jp   - m-element cell array of indices for p (indexed evaluatons)
%               Jy   - m-element cell array of indices for y (indexed evaluatons)
%               w    - m-element cell array of probability weights for noise variables
%               vreorder - reordering index vector to change order of operations

function ws = fprocess(p,parents,X,e,options)
order=[];
mergevec=[];
Ionly=false;
if exist('options','var') && ~isempty(options)
  if isfield(options,'order'),       order       = options.order;       end
  if isfield(options,'mergevec'),    mergevec    = options.mergevec;    end
  if isfield(options,'Ionly'),       Ionly       = options.Ionly;       end
end
m=length(p);
if nargin < 3, error('must pass 3 inputs'); end
% get information on random noise (e) variables
if nargin<4 || isempty(e), e = {};    end
if ~iscell(e), e={e};  end

% remove unused variables
[parents,X,e,dx,de] = removevariables(parents,X,e);

% reorder
vreorder = []; 
if ~isempty(order) && ~isequal(1:m,order(:)')  
  % create an index to reorder v
  pm = cellfun(@(x)size(x,1),p);
  vreorder = shuffleindex(pm,order,1);
  p = p(order);
  parents = parents(order);
end

% reorders e variables so first summed is listed first 
if isempty(e)
  eenter = []; eexit = []; evals = {};
else
  [parents,e] = ereorder(parents,e,mergevec);
  [eenter,eexit] = getenterexit(parents,e);
  evals = cellfun(@(x)x.values,e,'UniformOutput',false);
end
% merge cpts
if isempty(mergevec)
  mergevec = EVoptgroups(p,parents,X,e,options);
end
k = 0;
m = length(mergevec);
for i=1:m
  ind = k+1:k+mergevec(i);  % merge variables
  % list of noise variables that can be removed
  eremove = eenter>=ind(1) & eexit<=ind(end); 
  [p{i},parents{i}] = mergecpts(p(ind),parents(ind),X,e,eremove);
  k = ind(end);  
end
p = p(1:m);
parents = parents(1:m);
% remove unused variables
[parents,X,e] = removevariables(parents,X,e);
[eenter,eexit] = getenterexit(parents,e);
evals = cellfun(@(x)x.values,e,'UniformOutput',false);
ne = cellfun(@length,evals);

% get the weighting vectors for each stage
w = cell(1,m);
if ~isempty(e)
  eweights=cellfun(@(x)x.cpt,e,'UniformOutput',false);
  for i = m:-1:1
    ii = eexit == i;
    if any(ii)
      ii = find(ii);
      wi = eweights{ii(1)};
      for j=2:length(ii)
        wi = eweights{ii(j)}*wi';
        wi = wi(:);
      end
      w{i} = wi;
    end
  end
end

[Ip,Iy] = EVindices(parents,X,e,1);

if ~Ionly
  [Jp,Jy] = EVgetJ(Ip,Iy,eenter,eexit,ne);
else
  Jp=[];  Jy=[];
end

for i=1:m
  if Ip{i}(end) == numel(Ip{i}), Ip{i} = []; end
  if Iy{i}(end) == numel(Iy{i}), Iy{i} = []; end
end

% pack the workspace
ws.p=p;
ws.Ip=Ip;
ws.Iy=Iy;
ws.Jp=Jp;
ws.Jy=Jy;
ws.w=w;
ws.vreorder=vreorder;
return

% ereorder Finds the order of noise variables that has them sorted
% first by how they are summed out in a mergevec operation and then
% by how the remaining variables are summed out in the EV evaluation
%
% INPUTS
%   parents  : m-element cell array of conditioning variables (parents)
%   e        : cell array of rv structures (discrete or w/ discrete approximations)
%   mergevec : m-vector with the size of each group (sum(mergevec)=ds)
% OUTPUT
%   parents  : revised parents with new number of noise variables
%   e        : reordered noise variables
%   eorder   : new order for the noise variables
function [parents,e,eorder]=ereorder(parents,e,mergevec)
% obtain when each noise variable enters and exits
[eenter,eexit]=getenterexit(parents,e);
% handle 0 or 1 entering noise terms
ne = length(eenter);
eorder = [];
if ne<=1
  if ne==1, eorder = 1; end
  return
end
ds=length(parents);
m = length(mergevec);
notsummable = eexit>0;
k = 0;
for i = 1:m
  ind = k+1 : k+mergevec(i);
  summable = find(eenter>=ind(1) & eexit<=ind(end));
  [~,ii] = sort(eexit(summable));
  eorder = [eorder; summable(ii)]; %#ok<AGROW>
  notsummable(summable) = false;
  k = k + mergevec(i);
end
notsummable = find(notsummable); 
[~,ii] = sort(eexit(notsummable));
eorder = [eorder; notsummable(ii)]; 
for i=1:ds
  ii = parents{i}<0;
  parents{i}(ii) = -matchindices(-parents{i}(ii),eorder);
end
e=e(eorder);
  
    





function [parents,X,e,dx,de]=removevariables(parents,X,e)  
m=length(parents);
de = length(e);
if iscell(X)
  dx = length(X);
else
  dx = size(X,2);
end
% remove any unused variables
allparents = unique([parents{:}]);
alle = sort(-allparents(allparents<0));
dep = length(alle);
if dep < de
  for i=1:m
    ii=parents{i}<0;
    parents{i}(ii) = -matchindices(-parents{i}(ii),alle);
  end
  e=e(alle);
  de=dep;
end
allx = allparents(allparents>0);
dxp = length(allx);
if dxp < dx
  for i=1:m
    ii=parents{i}>0;
    parents{i}(ii) = matchindices(parents{i}(ii),allx);
  end
  if iscell(X)
    X=X(allx);
  else
    X = unique(X(:,allx),'rows');
  end
  dx=dxp;
end

% EVgetJ generate the J indices 
% USAGE
%    [Jp,Jy] = EVgetJ(Ip,Iy);
% INPUTS
%   Ip, Iy : index vectors for P and y at each stage relative to combined set
% OUTPUTS
%   Jp, Jy : index vectors for P and y at each stage relative to X
function [Jp,Jy] = EVgetJ(Ip,Iy,eenter,eexit,ne)
if nargin<3 || isempty(eexit)
  nonoise = true;
else
  nonoise = false;
end
m = length(Iy);
Jp = cell(1,m); 
Jy = cell(1,m);
nX = length(Iy{end});
if nonoise
  nem = 1;
else
  nem = prod(ne(eexit==m));
  nX = nX/nem;
end
J = (1:nX)';  
for k = m:-1:1
  kk = k>=eenter & k<=eexit;  
  nek = prod(ne(kk));
  if ~isempty(Ip{k}), 
    Jp{k}=reshape(Ip{k},[],nek);
    Jp{k} = Jp{k}(J,:); 
  end
  if ~isempty(Iy{k}), 
    Jy{k}=reshape(Iy{k},[],nek);
    Jy{k} = Jy{k}(J,:); 
    J = Jy{k}(:,1);
  end
  if Jp{k}(end) == numel(Jp{k}), Jp{k} = []; end
  if Jy{k}(end) == numel(Jy{k}), Jy{k} = []; end
end


% EVeval Evaluates an EV function using EVmergefunc
% USAGE
%   y=EVeval(p,Ip,Iy,Jp,Jy,vreorder,useI,alg,firstJ,v,Ie);
% INPUTS
%   p  : 1 x m cell array of conditional probability matrices
%                ith element has n(i) rows
%   Ip : 1 x m cell array of expansion indices for the p inputs
%   Iy : 1 x m cell array of expansion indices for the output at each stage
%   Jp : 1 x m cell array of expansion indices at each stage relative to X
%   vreorder : index vector to rearrange v (empty implied no rearrangement)
%   useI     : 0 to force J indexing at the first iteration
%   alg      : not used
%   firstJ   : first step using J indices
%   v        : prod(n)-vector
%   Ie       : index of extraction values relative to X - used to evaluate 
%                the expected value conditioned on selected values of X
% OUTPUT
%   y  : E[v]

% Note: empty index vector avoids unnecessary indexing by indicating that all
%       columns are used
function y=EVeval(p,w,Ip,Iy,Jp,Jy,vreorder,useI,alg,firstJ,v,Ie)
if nargin>11, 
  extract=true;   nIe=length(Ie); 
else
  extract=false;
end
if ~isempty(vreorder), v=v(vreorder); end
m=length(p);
ni=size(p{1},1);

if ~extract  
  if isempty(Ip{1})  % check if p{1} needs to be reordered
     y = reshape(v,[],ni) * p{1};
  else
     y = reshape(v,[],ni) * p{1}(:,Ip{1});
  end
  for i=2:m
    y=indexedmult(y,Iy{i},p{i},Ip{i},false,alg,w{i});
  end
else
  % p{1} is not handled correctly in the else clause so avoid it for now
  if ((useI && nIe*size(Jp{1},2)>=size(p{1},2)) || firstJ>1) && (m>1)
    if isempty(Ip{1})  % check if p{1} needs to be reordered 
      y = reshape(v,[],ni) * p{1}; 
    else
      y = reshape(v,[],ni) * p{1}(:,Ip{1}); 
    end
    expanded=false;
  else
    if isempty(Jp{1}),   pind=uint64(Ie);
    else                 pind=Jp{1}(Ie,:);
    end
    y = reshape(v,[],ni) * p{1}(:,pind);
    useI=false;
    expanded=true;
  end
  for i=2:m
    % determine if switchover to J indexing should occur (if it hasn't already)
    if useI && (nIe*size(Jp{i},2)<max(length(Iy{i}),size(y,2)) || i==m) && i>=firstJ
      %disp(['Using J indexing starting in iteration ' num2str(i)])
      useI=false;
    end
    if useI
      y=indexedmult(y,Iy{i},p{i},Ip{i},false,alg,w{i});  
      if i==m, y = y(Ie); end
    else
      if expanded, 
        yind=[];
      else
        if isempty(Jy{i}), yind=uint64(Ie);
        else               yind=Jy{i}(Ie,:);
        end
      end
      if isempty(Jp{i}),   pind=uint64(Ie);
      else                 pind=Jp{i}(Ie,:);
      end
      y=indexedmult(y,yind,p{i},pind,false,alg,w{i});
      expanded=true;
    end
  end
end
y=y(:);
