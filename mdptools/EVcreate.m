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
%   usebsxfun : forces the use of bsxfun rather than tprod
%   useI      : if 0 uses J indexing from the beginning
%
% Note: in general the workspace output is not needed but is provided here for  
% debugging purposes or to help curious users understand how this function works
function [EV,workspace]=EVcreate(p,X,parents,options)
if nargin==2
  ws=X;
  if ~isempty(ws.order)
    m=length(p);
    n=cellfun(@(x)size(x,1),p);
    p=p(ws.order);
    % create an index to reorder v
    vreorder=reshape(1:prod(n),n(m:-1:1));
    vreorder=permute(vreorder,m+1-ws.order(m:-1:1));
    vreorder=vreorder(:);
  end
  if exist('tprod','file') 
    EV=@(varargin) evalEV(p,ws.Ip,ws.Iy,ws.Jp,ws.Jy,vreorder,1,varargin{:}); 
  else 
    EV=@(varargin) evalEVb(p,ws.Ip,ws.Iy,ws.Jp,ws.Jy,vreorder,1,varargin{:}); 
  end
else
  m=length(p);
  order=[];
  usebsxfun=false;
  useI=true;
  if exist('options','var') && ~isempty(options)
    if isfield(options,'order'),     order     = options.order;     end
    if isfield(options,'usebsxfun'), usebsxfun = options.usebsxfun; end
    if isfield(options,'useI'),      useI   = options.useI;         end
  end
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
  else
    vreorder=[];
  end
  Ip=cell(1,m); Iy=cell(1,m); Jp=cell(1,m); Jy=cell(1,m);
  parentsall=unique([parents{:}]);
  parentsout=parents{1};
  [Xout,~,Jp{1}]=unique(X(:,parents{1}),'rows');
  if size(Xout,1)~=size(p{1},2)
    error('parents{1} is incompatible with p{1}')
  end
  for i=2:m
    parentscombined=union(parentsout,parents{i});
    Xcombined=unique(X(:,parentscombined),'rows');
    if isequal(parents{i},parentscombined) 
      Ip{i}=[]; 
      mi=size(Xcombined,1);
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
  if exist('tprod','file') && ~usebsxfun
    EV=@(varargin) evalEV(p,Ip,Iy,Jp,Jy,vreorder,useI,varargin{:}); 
  else 
    EV=@(varargin) evalEVb(p,Ip,Iy,Jp,Jy,vreorder,useI,varargin{:}); 
  end
  if nargout>1
    workspace.p=p;
    workspace.Ip=Ip;
    workspace.Iy=Iy;
    workspace.Jp=Jp;
    workspace.Jy=Jy;
    workspace.vreorder=vreorder;
  end
end


% Note: empty index vector avoid unnecessary indexing
function y=evalEVt(p,Ip,Iy,Jp,Jy,vreorder,useI,v,Ie)
if nargin>8, extract=true; nIe=length(Ie); else extract=false; end
if ~isempty(vreorder), v=v(vreorder); end
m=length(p);
ni=size(p{1},1);
if ~extract || (useI && nIe>=size(p{1},2))
  y = reshape(v,length(v)/ni,ni) * p{1}; 
else
  y = reshape(v,numel(v)/ni,ni) * p{1}(:,Jp{1}(Ie));
  %disp('Using J indexing starting in iteration 1')
  %   the following line seems like it should be faster than the one above
  %   but indexing is slow and expanding columns of y, which is 
  %   generally much bigger than p{1}, slows down this operation
  %y = reshape(v,length(v)/ni,ni) * p{1}; y = y(:,Jp{1}(Ie));
  useI=false;
end
for i=2:m
  % determine if switchover to J indexing should occur (if it hasn't already)
  if extract && useI && (nIe<=max(length(Iy{i}),size(y,3)) || i==m)
    %disp(['Using J indexing starting in iteration ' num2str(i)])
    useI=false;
    if isempty(Jy{i}),   y=y(:,Ie);
    else                 y=y(:,Jy{i}(Ie));
    end
  end
  ni=size(p{i},1);
  y=reshape(y,[size(y,1)/ni,ni,size(y,2)]);
  if useI
    if ~isempty(Iy{i}), y=y(:,:,Iy{i}); end
    if ~isempty(Ip{i}), y=tprod(y, [1 -1 2], p{i}(:,Ip{i}),[-1 2]);
    else                y=tprod(y, [1 -1 2], p{i},         [-1 2]);
    end
  else
    if ~isempty(Jp{i}), y=tprod(y, [1 -1 2], p{i}(:,Jp{i}(Ie)), [-1 2]);
    else                y=tprod(y, [1 -1 2], p{i}(:,Ie),        [-1 2]);
    end
  end
end
y=y(:);


% EVevalt Evaluates an EV function using tprod
% USAGE
%   y=EVevalt(p,Ip,Iy,Jp,Jy,vreorder,useI,v,Ie);
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

% Note: empty index vector avoid unnecessary indexing
function y=evalEV(p,Ip,Iy,Jp,Jy,vreorder,useI,v,Ie)
if nargin>8, extract=true; nIe=length(Ie); else extract=false; end
if ~isempty(vreorder), v=v(vreorder); end
m=length(p);
ni=size(p{1},1);
if ~extract || (useI && nIe>=size(p{1},2))
  y = reshape(v,length(v)/ni,ni) * p{1}; 
else
  y = reshape(v,numel(v)/ni,ni) * p{1}(:,Jp{1}(Ie));
  %disp('Using J indexing starting in iteration 1')
  %   the following line seems like it should be faster than the one above
  %   but indexing is slow and expanding columns of y, which is 
  %   generally much bigger than p{1}, slows down this operation
  %y = reshape(v,length(v)/ni,ni) * p{1}; y = y(:,Jp{1}(Ie));
  useI=false;
end
for i=2:m
  % determine if switchover to J indexing should occur (if it hasn't already)
  if extract && useI && (nIe<=max(length(Iy{i}),size(y,3)) || i==m)
    %disp(['Using J indexing starting in iteration ' num2str(i)])
    useI=false;
  end
  if useI
    y=EVmergefunc(y,Iy{i},p{i},Ip{i});
  else
    if isempty(Jy{i}),   yind=uint64(Ie);
    else                 yind=Jy{i}(Ie);
    end
    if isempty(Jp{i}),   pind=uint64(Ie);
    else                 pind=Jp{i}(Ie);
    end
    y=EVmergefunc(y,yind,p{i},pind);
  end
end
y=y(:);


function y=evalEVb(p,Ip,Iy,Jp,Jy,vreorder,useI,v,Ie)
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

