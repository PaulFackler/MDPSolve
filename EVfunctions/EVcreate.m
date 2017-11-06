% EVcreate Creates an EV function from conditional probabilities
% USAGE
%   EV=EVcreate(p,X,parents,options);
% INPUTS
%   p       : m-element cell array of conditional probability matrices
%   X       : matrix of state/action combinations
%   parents : m-element cell array of conditioning variables (parents)
%   options : structure variable (fields described below)
% OUTPUTS
%   EV        :  a function handle for an EV function
%   workspace : a structure variable containing the workspace of the EV function
% parents{i} contains the columns of X that condition the ith variable
% Thus Xi=unique(X(:,parents{i}),'rows') contains the unique elements of the
%   conditioning variables that condition output variable i.
% It must be the case that the number of rows in Xi matches the number 
%   of columns in p{i}: size(Xi,1)=size(p{i},2)
%
% Options fields:
%   order     : m-vector if alternative order of operation
%   usebsxfun : forces the use of bsxfun rather than EVmergefunc
%   useI      : if 0 uses J indexing from the beginning
%   vreorder  : reordering vector for v (ignored if order option is defined)
%   mergevec  : vector providing merge information with sum(mergevec)=m
%                 An EV function is created with length(mergevec) components
%                 The default is ones(1,m), i.e., one component for each variable
%                 mergevec=m results in a single transition matrix
%   getoptord : attempts to find optimal combination method
%   expandall : expands all of the CPTS to have nx columns
%   nx        : sizes of the X variables [used to determine optimal grouping]
%
% Note: in general the workspace output is not needed but is provided here for  
% debugging purposes or to help curious users understand how this function works
%
% If EVcreate is called multiple times using different values of p but with
%   no change to X, parents or options, then call it the first time using
%     [EV,workspace]=EVcreate(p,X,parents,options);
% and on subsequent calls use
%     EV=EVcreate(p,workspace);
function [EV,workspace]=EVcreate(p,X,parents,options)
if nargin==2
  ws=X;
  % check if v needs to be reordered
  if ~isempty(ws.order) && any(diff(ws.order)~=1)
    m=length(p);
    n=cellfun(@(x)size(x,1),p);
    p=p(ws.order);
    % create an index to reorder v
    vreorder=reshape(1:prod(n),n(m:-1:1));
    vreorder=permute(vreorder,m+1-ws.order(m:-1:1));
    vreorder=vreorder(:);
  else
    vreorder = [];
  end
  if exist('EVmergefunc','file') 
    EV=@(varargin) EVeval(p,ws.Ip,ws.Iy,ws.Jp,ws.Jy,vreorder,1,varargin{:}); 
  else 
    EV=@(varargin) EVevalb(p,ws.Ip,ws.Iy,ws.Jp,ws.Jy,vreorder,1,varargin{:}); 
  end
  workspace=ws;
else
  m=length(p);
  order=[];
  usebsxfun=false;
  useI=true;
  vreorder=[];
  mergevec=[];
  getoptord=false;
  expandall=false;
  mm=[];
  nx=[];
  if exist('options','var') && ~isempty(options)
    if isfield(options,'order'),     order      = options.order;     end
    if isfield(options,'usebsxfun'), usebsxfun  = options.usebsxfun; end
    if isfield(options,'useI'),      useI       = options.useI;      end
    if isfield(options,'vreorder'),  vreorder   = options.vreorder;  end
    if isfield(options,'mergevec'),  mergevec   = options.mergevec;  end
    if isfield(options,'mm'),        mm         = options.mm;        end
    if isfield(options,'getoptord'), getoptord  = options.getoptord; end
    if isfield(options,'expandall'), expandall  = options.expandall; end
    if isfield(options,'nx'),        nx         = options.nx;        end
  end
  if getoptord
    [p,parents,order,vreorder]=EVpreprocess(p,X,parents);
    m=length(p);
  else
    % check if mergevec is correctly specified
    if ~isempty(mergevec) && sum(mergevec)~=m,
      error('mergevec is incorrectly specified')
    end
    if length(mergevec)==m, mergevec=[]; end  % no preprocessing needed

    if ~isempty(order)
      if length(order)~=m
        error(['options.order must be an ' num2str(m) '-vector']);
      end
      n=cellfun(@(x)size(x,1),p);
      p=p(order);
      parents=parents(order);
      % create an index to reorder v
      vreorder=reshape(1:prod(n),n(m:-1:1));
      vreorder=permute(vreorder,m+1-order(m:-1:1));
      vreorder=vreorder(:);
    end
    % get optimal grouping
    if isempty(mergevec) && ~isempty(mm)
      pp=cellfun(@(x)size(x,1),p);
      pp=fliplr(cumprod(fliplr(pp)));
      mergevec=optmergeorder(pp,mm,options);
    end

    % combine CPTs as needed
    if ~isempty(mergevec)
      k=0;
      m=length(mergevec);
      for i=1:m
        [p{i},parents{i}]=mergecpts(p(k+1:k+mergevec(i)),parents(k+1:k+mergevec(i)),X);
        k=k+mergevec(i);
      end
      % check if full transition matrix is used
      p=p(1:m);
      parents=parents(1:m);
      if m==1
        EV=@(varargin) EVuseP(p{1},varargin{:});
        return
      end
    end
  end
  
  % create index vectors
  Ip=cell(1,m); Iy=cell(1,m); Jp=cell(1,m); Jy=cell(1,m);
  mm=zeros(1,m);
  parentsall=unique([parents{:}]);
  parentsout=parents{1};
  [Xout,~,Jp{1}]=unique(X(:,parents{1}),'rows');
  if size(Xout,1)~=size(p{1},2)
    error('parents{1} is incompatible with p{1}')
  end
  mm(1)=size(Xout,1);
  for i=2:m
    parentscombined=union(parentsout,parents{i});
    if ~isequal(parentsout,parentscombined)
      Xcombined=unique(X(:,parentscombined),'rows');
    end
    mm(i)=size(Xcombined,1);
    if isequal(parents{i},parentscombined) 
      Ip{i}=[]; 
      mi=mm(i);
    else
      ii=ismember(parentscombined,parents{i});
      [~,~,Ip{i}]=unique(Xcombined(:,ii),'rows');
      mi=max(Ip{i});
    end
    if mi~=size(p{i},2)
      error(['parents{' num2str(i) '} is incompatible with p{' num2str(i) '}'])
    end
    if length(parents{i})==length(parentsall)
      Jp{i}=[];
    else
      [~,~,Jp{i}]=unique(X(:,parents{i}),'rows');
    end
    if isequal(parentsout,parentscombined) 
      Iy{i}=[];
    else
      ii=ismember(parentscombined,parentsout);
      [~,~,Iy{i}]=unique(Xcombined(:,ii),'rows');
    end
    if length(parentsout)==length(parentsall)
      Jy{i}=[];
    else
      [~,~,Jy{i}]=unique(X(:,parentsout),'rows');
    end
    parentsout=parentscombined;
  end
  clear X Xin Xout Xcombined parentsin parentsout parentscombined i m
  if exist('EVmergefunc','file') && ~usebsxfun
    EV=@(varargin) EVeval(p,Ip,Iy,Jp,Jy,vreorder,useI,varargin{:}); 
  else 
    EV=@(varargin) EVevalb(p,Ip,Iy,Jp,Jy,vreorder,useI,varargin{:}); 
  end
  if nargout>1
    workspace.p=p;
    workspace.Ip=Ip;
    workspace.Iy=Iy;
    workspace.Jp=Jp;
    workspace.Jy=Jy;
    workspace.vreorder=vreorder;
    workspace.mm=mm;
    workspace.mergevec=mergevec;
  end
end

% special case when the full transition matrix is used
function y=EVuseP(P,v,Ix)
if nargin==2
  y=P'*v;
else
  y=P(:,Ix)'*v;
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
function y=EVeval(p,Ip,Iy,Jp,Jy,vreorder,useI,v,Ie)
if nargin>8, extract=true; nIe=length(Ie); else extract=false; end
if ~isempty(vreorder), v=v(vreorder); end
m=length(p);
ni=size(p{1},1);
if ~extract
  y = reshape(v,length(v)/ni,ni) * p{1}; 
  for i=2:m
    y=EVmergefunc(y,Iy{i},p{i},Ip{i});
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
      y=EVmergefunc(y,Iy{i},p{i},Ip{i});
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
      y=EVmergefunc(y,yind,p{i},pind);
      expanded=true;
    end
  end
end

y=y(:);


function y=EVevalb(p,Ip,Iy,Jp,Jy,vreorder,useI,v,Ie)
if nargin>8, extract=true; nIe=length(Ie); else extract=false; end
if ~isempty(vreorder), v=v(vreorder); end
m=length(p);
ni=size(p{1},1);
if ~extract || (useI && nIe>=size(p{1},2))
  y = reshape(v,length(v)/ni,ni) * p{1}; 
else
  y = reshape(v,numel(v)/ni,ni) * p{1}(:,Jp{1}(Ie));
  useI=false;
end
for i=2:m
  % determine if switchover to J indexing should occur (if it hasn't already)
  if extract && useI && (nIe<=max(length(Iy{i}),size(y,3)) || i==m)
    %disp(['Using J indexing starting in iteration ' num2str(i)])
    useI=false;
    nIe=length(Ie);
    if isempty(Jy{i}),   y=y(:,Ie);
    else                 y=y(:,Jy{i}(Ie));
    end
  end
  ni=size(p{i},1);
  y=reshape(y,[size(y,1)/ni,ni,size(y,2)]);
  if useI
    if ~isempty(Iy{i}), y=y(:,:,Iy{i});  end
    nIpi=length(Ip{i});
    if nIpi>0, y=bsxfun(@times,y,reshape(p{i}(:,Ip{i}),[1 ni nIpi]));
    else       y=bsxfun(@times,y,reshape(p{i},         [1 size(p{i})]));
    end
  else
    if ~isempty(Jp{i}), y=bsxfun(@times,y,reshape(p{i}(:,Jp{i}(Ie)),[1 ni nIe]));
    else                y=bsxfun(@times,y,reshape(p{i}(:,Ie),       [1 ni nIe]));
    end
  end
  y=squeeze(sum(y,2));
end
y=y(:);


function z=merge(x,xind,y,yind)
p=size(y,1);
n=size(x,1)/p;
if isempty(xind)
  m=size(x,2);
else
  m=length(xind);
end
x=reshape(x,[n,p,size(x,2)]);
z=zeros(n,m);
if isempty(xind)
  if isempty(yind)
    for i=1:m
      z(:,i)=x(:,:,i)*y(:,i);
    end
  else
    for i=1:m
      z(:,i)=x(:,:,i)*y(:,yind(i));
    end
  end
else
  if isempty(yind)
    for i=1:m
      z(:,i)=x(:,:,xind(i))*y(:,i);
    end
  else
    for i=1:m
      z(:,i)=x(:,:,xind(i))*y(:,yind(i));
    end
  end
end

