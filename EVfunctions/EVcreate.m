% EVcreate Creates an EV function from conditional probabilities
% USAGE
%   EV=EVcreate(p,X,parents,options);
% INPUTS
%   p       : m-element cell array of conditional probability matrices
%   X       : matrix of state/action combinations
%   parents : m-element cell array of conditioning variables (parents)
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
%   order       : m-vector if alternative order of operation
%   usebsxfun   : forces the use of bsxfun rather than EVmergefunc
%   useI        : if 0 uses J indexing from the beginning
%   vreorder    : reordering vector for v (ignored if order option is defined)
%   mergevec    : vector providing merge information with sum(mergevec)=m
%                   An EV function is created with length(mergevec) components
%                   The default is ones(1,m), i.e., one component for each variable
%                   mergevec=m results in a single transition matrix
%   getoptord   : attempts to find optimal combination method
%   getoptgroup : determines the optimal grouping for a given order
%   penalty     : penalty for more factors
%   nx          : sizes of the X variables (used to determine optimal grouping)
%   expandall   : expands all of the CPTS to have nx columns
%
% Note: in general the workspace output is not needed but is provided here for  
% debugging purposes or to help curious users understand how this function works
%
% If EVcreate is called multiple times using different values of p but with
%   no change to X, parents or options, then call it the first time using
%     [EV,workspace]=EVcreate(p,X,parents,options);
% and on subsequent calls use
%     EV=EVcreate(p,workspace);
%
% EVcreate can also be called as EV=EVcreate(P) or as EV=EVcreate({EV0,EV1})
% to facilitate defining an EV functions that does both full and indexed evaluations.
function [EV,ws]=EVcreate(p,X,parents,options)
if nargin==1 
  ws=[];
  if isnumeric(p)
    EV=@(varargin) EVuseP(p,varargin{:});
    return
  elseif iscell(p) && length(p)==2
    EV=@(varargin) EVcombine(p{1},p{2},varargin{:});
    return
  end
end
order=[];
usebsxfun=false;
useI=true;
vreorder=[];
getoptord=false;
mergevec=[];
pm=[];
nx=[];
expandall=false;
getoptgroup=true;
ws=[];
if exist('options','var') && ~isempty(options)
  if isfield(options,'order'),       order       = options.order;       end
  if isfield(options,'usebsxfun'),   usebsxfun   = options.usebsxfun;   end
  if isfield(options,'useI'),        useI        = options.useI;        end
  if isfield(options,'getoptord'),   getoptord   = options.getoptord;   end
  if isfield(options,'vreorder'),    vreorder    = options.vreorder;    end
  if isfield(options,'mergevec'),    mergevec    = options.mergevec;    end
  if isfield(options,'getoptgroup'), getoptgroup = options.getoptgroup; end
  if isfield(options,'pm'),          pm          = options.pm;          end
  if isfield(options,'nx'),          nx          = options.nx;          end
  if isfield(options,'expandall'),   expandall   = options.expandall;   end
  if isfield(options,'workspace'),   ws          = options.workspace;   end
else
  options=[];
end

m=length(parents);
% check if mergevec is correctly specified
if ~isempty(mergevec) && sum(mergevec)~=m,
  error('mergevec is incorrectly specified - must sum to m')
end
%if length(mergevec)==m, mergevec=[]; end  % no preprocessing needed


% check if order is correctly specified
if ~isempty(order)
  if length(order)~=m
    error(['options.order must be an ' num2str(m) '-vector']);
  end
  if ~isequal(1:m,sort(order(:))')
    error(['options.order must be a permutation of 1:' num2str(m)]);
  end
  if isequal(1:m,order(:)')
    order=[];  % no reordering needed
  end
end

if ~isempty(ws) 
  [Ip,Iy,Jp,Jy,Xp,yn,vreorder,mergevec]=unpackworkspace(ws);
end 

% get CPT sizes
if ~isempty(p)
  [pm,pn]=cellfun(@size,p);
else
  pn=[];
end

if getoptord
  [order,mergevec]=EVoptorder(p,X,parents,nx,options); 
  getoptgroup=false;
end
% reorder
if ~isempty(order)
  % create an index to reorder v
  vreorder=shuffleindex(pm,order,1);
  p=p(order);
  parents=parents(order);
  if isempty(mergevec) && ~getoptgroup
    if ~isempty(pn)
      pn=pn(order);
      [Ip,Iy,Jp,Jy,Xp,yn]=EVindices(X,parents,pn,options); % performs checks on p
    else
      [Ip,Iy,Jp,Jy,Xp,yn]=EVindices(X,parents,[],options);           % no checks
    end
  end
end
% get optimal grouping
if isempty(mergevec) && getoptgroup
  pm=fliplr(cumprod(fliplr(pm)));
  if ~exist('yn','var') || isempty(yn)
    yn=EVgetyn(X,parents,options);
  end
  mergevec=EVoptgroups(pm,yn,options);
end

% combine CPTs as needed
if ~isempty(mergevec) && length(mergevec)<length(parents)
  k = 0;
  m = length(mergevec);
  for i=1:m
    ind = k+1:k+mergevec(i);
    %tic
    %[p2{i},parents2{i}] = mergecpts2(p(ind),X,parents(ind));
    %tt=toc;
    %tic
    [p{i},parents{i}] = mergecpts(p(ind),X,parents(ind),options);
    %disp(toc/tt)
    %if ~isequal(p{i},p2{i})
    %  error('oops')
    %end
    k = ind(end);
  end
  p = p(1:m);
  parents = parents(1:m);
  [Ip,Iy,Jp,Jy,Xp,yn] = EVindices(X,parents,[],options); 
end

% check if full transition matrix is used
if m==1
  p=p{1};
  EV=@(varargin) EVuseP(p,varargin{:});
  return
end

% if indices are not yet obtained then get them now
if ~exist('Ip','var') 
  if ~isempty(p)
    pn=cellfun(@(x)size(x,2),p);
    [Ip,Iy,Jp,Jy,Xp,yn]=EVindices(X,parents,pn,options); % performs checks on p
  else
    [Ip,Iy,Jp,Jy,Xp,yn]=EVindices(X,parents,[],options);    % no checks
  end
end

clear X expandall getoptord i m nx order pm
EV=@(varargin) EVeval(p,Ip,Iy,Jp,Jy,vreorder,useI,usebsxfun,varargin{:}); 
if nargout>1
  ws=packworkspace(Ip,Iy,Jp,Jy,Xp,yn,vreorder,mergevec);
end


% packworkspace Creates a workspace structure
function workspace=packworkspace(Ip,Iy,Jp,Jy,Xp,yn,vreorder,mergevec)
  workspace.Ip=Ip;
  workspace.Iy=Iy;
  workspace.Jp=Jp;
  workspace.Jy=Jy;
  workspace.Xp=Xp;
  workspace.yn=yn;
  workspace.vreorder=vreorder;
  workspace.mergevec=mergevec;
  
% unpackworkspace Extracts the elements of a workspace structure
function [Ip,Iy,Jp,Jy,Xp,yn,vreorder,mergevec]=unpackworkspace(workspace)
  Ip=workspace.Ip;
  Iy=workspace.Iy;
  Jp=workspace.Jp;
  Jy=workspace.Jy;
  Xp=workspace.Xp;
  yn=workspace.yn;
  vreorder=workspace.vreorder;
  mergevec=workspace.mergevec;
  

% special case when the full transition matrix is used
function y=EVuseP(P,v,Ix)
if nargin==2
  y=P'*v;
else
  y=P(:,Ix)'*v;
end

function y=EVcombine(EV0,EV1,v,Ix)
if nargin==3
  y=EV0(v);
else
  y=EV1(v,Ix);
end

% EVeval Evaluates an EV function using EVmergefunc
% USAGE
%   y=EVeval(p,Ip,Iy,Jp,Jy,vreorder,useI,v,Ie);
% INPUTS
%   p  : 1 x m cell array of conditional probability matrices
%                ith element has n(i) rows
%   Ip : 1 x m cell array of expansion indices for the p inputs
%   Iy : 1 x m cell array of expansion indices for the output at each stage
%   Jp : 1 x m cell array of expansion indices at each stage relative to X
%   vreorder : index vector to rearrange v (empty implied no rearrangement)
%   useI     : 0 to force J indexing at the first iteration
%   v        : prod(n)-vector
%   Ie       : index of extraction values relative to X - used to evaluate 
%                the expected value conditioned on selected values of X
% OUTPUT
%   y  : E[v]

% Note: empty index vector avoids unnecessary indexing by indicating that all
%       columns are used
function y=EVeval(p,Ip,Iy,Jp,Jy,vreorder,useI,usebsxfun,v,Ie)
if nargin>9, extract=true; nIe=length(Ie); else extract=false; end
if ~isempty(vreorder), v=v(vreorder); end
m=length(p);
ni=size(p{1},1);
if ~extract
  y = reshape(v,length(v)/ni,ni) * p{1}; 
  for i=2:m
    %y=EVmergefunc(y,Iy{i},p{i},Ip{i});
    y=indexedmult(y,Iy{i},p{i},Ip{i},false,usebsxfun);
  end
else
  if useI && nIe>=size(p{1},2)
    y = reshape(v,length(v)/ni,ni) * p{1}; 
    expanded=false;
  else
    if isempty(Jp{1}),   pind=uint64(Ie);
    else                 pind=Jp{1}(Ie);
    end
    y = reshape(v,numel(v)/ni,ni) * p{1}(:,pind);
    useI=false;
    expanded=true;
  end
  for i=2:m
    % determine if switchover to J indexing should occur (if it hasn't already)
    if useI && (nIe<=max(length(Iy{i}),size(y,3)) || i==m)
      %disp(['Using J indexing starting in iteration ' num2str(i)])
      useI=false;
    end
    if useI
      y=indexedmult(y,Iy{i},p{i},Ip{i},false,usebsxfun);
    else
      if expanded, yind=[];
      else
        if isempty(Jy{i}), yind=uint64(Ie);
        else               yind=Jy{i}(Ie);
        end
      end
      if isempty(Jp{i}),   pind=uint64(Ie);
      else                 pind=Jp{i}(Ie);
      end
      y=indexedmult(y,yind,p{i},pind,false,usebsxfun);
      expanded=true;
    end
  end
end
y=y(:);

% EVgetyn gets the sizes of the cumulative combined parent variables
% USAGE
%   yn=EVgetyn(X,parents,pn);
% INPUTS
%   X       : matrix of state/action combinations
%   parents : m-element cell array of conditioning variables (parents)
%   pn      : m-vector of column sizes of the CPTs (optional to check compatibility)
% OUTPUTS
%   yn          : m-vector of # of columns in the output at each evaluation step 
%
% parents{i} contains the columns of X that condition the ith variable
% Thus Xpi=unique(X(:,parents{i}),'rows') contains the unique elements of the
%   conditioning variables that condition output variable i.
% It must be the case that the number of rows in Xi matches the number 
%   of columns in p{i}: size(Xpi,1)=np(i)
  
function yn=EVgetyn(X,parents,options)
  m=length(parents);
  parentscombined=cell(1,m);
  parentscombined{1}=parents{1};
  for i=2:m
    parentscombined{i}=union(parentscombined{i-1},parents{i});
  end
  yn=zeros(1,m);
  Xi=X;
  for i=m:-1:2
    yn(i)=size(Xi,1);
    if ~isequal(parentscombined{i-1},parentscombined{i})
      ind=matchindices(parentscombined{i-1},parentscombined{i});
      [~,Xi]=getI(Xi,ind,options);
    end
  end
  yn(1)=size(Xi,1);
  
    
% matchindices Finds the values of ind1 that match those of ind0
function ind=matchindices(ind0,ind1)
  [~,ind]=ismember(ind0,ind1);
  ind=ind(ind>0);



% THIS DOES NOT WORK CORRECTLY
% combines index vectors for group merges  
function [Ip,Iy,Jp,Jy,pn]=combineindices(mergevec,Ip,Iy,Jp,Jy,pn)
cmv=cumsum(mergevec);
m=length(mergevec);
% need to handle missing index vectors (which indicate that parents and 
%   parentscombined are identical)
Jp{1}=Jp{1};
Jy{1}=ones(size(Jp{1}));
for i=1:m
  if i==1, lb=1;
  else     lb=cmv(i-1)+1;
  end
  Jp{i}=Jp{lb};
  Jy{i}=Jy{lb};
  Ipi=Ip{cmv(i)};
  Iyi=Iy{cmv(i)};
  for j=cmv(i)-1:-1:lb
    Ipi=Ip{j}(Ipi);
    Iyi=Iy{j}(Iyi);
  end
  Ip{i}=Ipi;
  Iy{i}=Iyi;
end
Ip=Ip(1:m);
Jp=Jp(1:m);
Iy=Iy(1:m);
Jy=Jy(1:m);
if nargin>=6, pn=pn(cmv); end
