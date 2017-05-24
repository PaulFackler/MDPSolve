% rvdef Defines a random number structure for a specified distribution
% USAGE 
%   rv=rvdef(type,parameters,values,weights);
% INPUTS
%   type        : distribution type (see availability below)
%   parameters  : parameter values or 
%                    a function that returns conditional parameter values
%   values      : a column n-vector of values or a scalar value for n (values
%                    and weights provided by the procedure when appropriate)
%   weights     : an n-row matrix of probability values (columns sum to 1)
% OUTPUT
%    rv         : a random variable (rv) structure with fields:
%                    type       : string indicating the distribution family
%                    parameters : parameter values
%                    values     : n-vector of values
%                    cpt        : conditional probability table
%                    size       : the value n
%                    order      : how CPT is organized
%                 (not all distributions use every field)  
%
% Available distributions:
%   'f'        deterministic         handle to function of parents
%   'c'        continuous variable   [lower bound;upper bound]
%   'd'        discrete              probability weights
%   'logit'    logit                 logit parameters
%   'i'        integer               [a;b] integer values on [a b] w/ equal weights
%   'u'        uniform               [a;b] range (with equal weights)
%   'ug'       uniform               [a;b] same as 'u' but uses Gaussian quadrature nodes & weights 
%   'n'        normal (Gaussian)     [mean;standard deviation] Gaussian quadrature nodes & weights 
%   'ne'       normal (Gaussian)     [mean;standard deviation] evenly spaced values
%   'g'        Gamma                 [a;b]=[shape;scale]       mean=a*b
%   'b'        Beta                  [a;b]                     mean=a/(a+b)
%   'burr12'   Burr-12               [a;b]
%   'k'        Kumaraswamy           [a;b]
%   'lin'      (1-c)+2cx on [0,1]    [c]                       mean=0.5+c/6 
%   'tri'      triangular            lower, upper, mode
%
% The parameters field can be a constant vector or can be a function handle
% This function should accept values of the parents and return values
% of the parameters (with rows corresponding to parameters and columns
% to realizations of the parent variables).
%
% Used by dsim and getprMC. Also can be used with rvgen to generate random
% values.
%
% for discrete distributions
%   order   : 'lc', 'lr', 'rc', 'rr'
%               'l' indicates lexicographic ordering of parent variables
%               'r' indicates reverse lexicographic ordering
%               'c' indicates column stochastic (children associated with rows)
%               'r' indicates row stochastic form
%               default: 'lc'

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2014, Paul L. Fackler (paul_fackler@ncsu.edu)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without  
% modification, are permitted provided that the following conditions are met:
% 
%    * Redistributions of source code must retain the above copyright notice, 
%        this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice, 
%        this list of conditions and the following disclaimer in the 
%        documentation and/or other materials provided with the distribution.
%    * Neither the name of the North Carolina State University nor of Paul L. 
%        Fackler may be used to endorse or promote products derived from this 
%        software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% For more information, see the Open Source Initiative OSI site:
%   http://www.opensource.org/licenses/bsd-license.php

function rv=rvdef(type,varargin)
switch type
  case 'f'
    rv=fdef(varargin{:});
  case 'v'
    rv=vdef(varargin{:});
  case 'd'
    rv=ddef(varargin{:});
  case 'n'
    rv=namedcrv('n',2,[-inf,inf],varargin{:});
  case 'ne'
    rv=namedcrv('ne',2,[-inf,inf],varargin{:});
  case {'u','ug'}
    if isempty(varargin) || isempty(varargin{1}), varargin{1}=[0;1]; end
    rv=namedcrv(type,2,[],varargin{:});
  case 'i'  % integer valued with even weights
    if length(varargin)~=1
      error('integer type rvs have only one input');
    end
    rv=rvint(varargin{1});
  case 'g'
    rv=namedcrv('g',2,[0,inf],varargin{:});
  case 'b'
    rv=namedcrv('b',2,[0,1],varargin{:});
  case 'burr12'
    rv=namedcrv('burr12',2,[0,inf],varargin{:});
  case 'k'
    rv=namedcrv('k',2,[0,1],varargin{:});
  case 'lin'
    rv=namedcrv('lin',1,[0,1],varargin{:});
  case 'tri'
    rv=namedcrv('tri',3,[],varargin{:});
  case 'logit'
    rv=logitdef(varargin{:});
  case 'bin'
    rv=bindef(varargin{:});
  case 'hypgeo'
    rv=hypgeodef(varargin{:});
  otherwise
    fname=[type 'def'];
    if exist(fname,'file')
      rv=feval(fname,varargin{:});
    else
      error(['distribution type ' type ' not available'])
    end
end



%%%%%%%%%%%%%%%%% UTILITY FUNCTIONS
% check if parameters matrix has the right number of values
function checksize(parameters,nr)
  if size(parameters,1)~=nr
    error(['parameters should have ' num2str(nr) ' rows'])
  end
  
% allocates inputs to outputs and fills no-passed inputs with emptys
function [varargout]=getinputs(mininputs,varargin)
nin=length(varargin);
if nin<mininputs
  error('not enough inputs passed')
end
if nargout<nin
  error('too many inputs passed')
end
varargout=cell(1,nargout);
for i=1:nin
  varargout{i}=varargin{i};
end
for i=nin+1:nargout
  varargout{i}=[];
end
    

% deterministic function
function rv=fdef(varargin)
[fh,values]=getinputs(1,varargin{:});
if isa(fh,'function_handle')    
  rv=struct('type','f','valfunc',fh);  
else
  error('for type f variables the second arguement should be a function handle')
end
if ~isempty(values)
  if size(values,2)>1
    error('for type f variables values should be a colummn vector')
  end
  rv.values=values;
end

% variable with values only
function  rv=vdef(varargin)
[lim,values]=getinputs(2,varargin{:});
if isempty(lim)
  lim=[-inf;inf];
else
  if length(lim)~=2
    error('lim should be a 2-element vector')
  elseif lim(1)>lim(2)
    error('lim(1) must be no greater than lim(2)')
  end
end
if isnumeric(values)
  if size(values,2)~=1
    error('values should be a column vector')
  else
    n=length(values);
    if n==0
      error('values should contain at least one element')
    end
  end
else
  error('values should be a column vector')
end
rv=struct('type','v','values',values,'size',n);
if lim(1)>-inf; rv.lb=lim(1); end;
if lim(2)< inf; rv.ub=lim(2); end

%%%%%%%%%%%%%%%%% DISCRETE RVs
% discrete rv defined by a CPT
function  rv=ddef(varargin)
switch nargin
  case 0
    error('must pass at least one input for d type variable')
  case 1
    params=varargin{1};
    values=[];
    order='cl';
  case 2
    params=varargin{1};
    values=varargin{2};
    order='cl';
  case 3
    params=varargin{1};
    values=varargin{2};
    order=varargin{3};
  otherwise      
    error('too many inputs for type d')
end
if isnumeric(params) 
  cpt=params;
  params=[];
  if any(order=='c')
    n=size(cpt,1);
    if max(abs(sum(cpt,1)-1))>1e-14
      warning('columns of parameter matrix do not sum to 1')
    end
  else
    n=size(cpt,2);
    if max(abs(sum(cpt,2)-1))>1e-14
      warning('rows of parameter matrix do not sum to 1')
    end
  end
  if any(any(cpt<0))
    warning('parameter matrix contains negative elements')
  end
else
  cpt=[];
  n=[];
end
if ~isempty(values) && isempty(n)
  n=length(values);
elseif ~isempty(values) && length(values)~=n
  error('values vector and parameter matrix are incompatible')  
end
if isempty(order)
  order='lc'; 
end
switch order
  case 'l'
    order='lc';
  case 'r'
    order='rc';
  case 'cl'
    order='lc';
  case 'cr'
    order='rc';
  case 'rl'
    order='lr';
  case []
    order='lc';
  case {'lc','lr','rc','rr'}
  otherwise
    error('order indicator must be ''lc'',''lr'',''rc'' or ''rr''')  
end
rv=struct('type','d','parameters',params,'values',values,'cpt',cpt,'order',order,'size',n);

function rv=logitdef(varargin)
[params,values]=getinputs(1,varargin{:});
n=size(params,1);
if isempty(values)
  values={1:n}';
else
  if size(values,1)~=n
    error('values vector not correctly specified for logit type variable')
  end
end
rv=struct('type','logit','parameters',params,'values',values,'size',n);

function rv=bindef(varargin)
[params,values]=getinputs(1,varargin{:});
if length(values)==1
  n=values;
  values=(0:n)';
else
  error('values field should be a scalar equal to the number of trials')
end
if iscell(params)
  if length(params)==1
    N=n;
  elseif length(params)==2
    N=params{2};
  else
    error('parameters incorrectly specified')
  end
  p=params{1};
elseif isnumeric(params)
  if length(params)==1
    N=n;
  elseif length(params)==2
    N=params(2);
  else
    error('parameters incorrectly specified')
  end
  p=params(1);
elseif isa(params,'function_handle')
  error('parameters cannot be defined by a function')
end
x=crectgrid(p,N,values);
P=binprob(x{:});
P=reshape(P,n+1,length(params)*length(N));
P(1,all(P==0,1))=1;
rv=struct('type','d','parameters',params,'values',values,'cpt',P,'size',n+1);

function rv=hypgeodef(varargin)
[params,values]=getinputs(1,varargin{:});
if length(values)==1
  n=values;
  values=(0:n)';
else
  error('values field should be a scalar equal to the number of trials')
end
if length(params)==1
  N=params;
else
  error('values field should be a scalar equal to the number of trials')
end
x=crectgrid((0:N)',(0:N)',values,values);
P=hypergeometric(x{:});
P=reshape(P,n+1,(N+1)^2*(n+1));
P(1,all(P==0,1))=1;
matchfunc=@(NN,KK,nn) NN+(KK+nn*(N+1))*(N+1)+1;
rv=struct('type','d','parameters',params,'values',values,'cpt',P,'size',N+1,'matchfunc',matchfunc);

% integer valued with uniform weights  
function rv=rvint(parameters) 
  if length(parameters)==1, parameters=[1;parameters]; end
  if parameters(2)<parameters(1)
    error('parameters(2)<parameters(1) for integer rvs is not allowed')
  end
  x=(parameters(1):parameters(2))';
  n=length(x);
  w=ones(n,1)/n;  
  rv=struct('type','i','parameters',parameters,'values',x,'cpt',w,'size',n);


%%%%%%%%%%%%%%%%% CONTINUOUS RVs
% named continuous distribution with discretization support
function rv=namedcrv(type,nparam,lim,varargin)
[params,values]=getinputs(1,varargin{:});
if isempty(values)
  n=0;
  values={[],[]};
elseif iscell(values)
  n=length(values{1});
  %  need better error checking !!!!!!!!!!!!!!!!!!
  if any(cellfun(@(x) any(size(x)~=[n 1]),values))
     error('values not correctly specified - should be a cell array with 2 column vectors')
  end
else
  n=values;
  if isnumeric(n) && all(size(n)==1) 
    values=cell(1,2);
    if n>0, [values{1},values{2}]=feval([type 'nw'],n,params);
    else    values={[],[]};
    end
  else
    error('second input should be a positive integer or a set of nodes and weights')
  end
end
checksize(params,nparam)
rv=struct('type',type,'parameters',params,'values',values{1},'cpt',values{2},'size',n);
if ~isempty(lim)
  if lim(1)>-inf, rv.lb=lim(1); end
  if lim(2)< inf, rv.ub=lim(2); end
end

function [x,w]=nnw(n,parameters) %#ok<*DEFNU>
  [x,w]=qnwnorm(n,0,1);
  x=bsxfun(@plus,x*parameters(2,:),parameters(1,:));
  
function [x,w]=nenw(n,parameters)
  [x,w]=qnwnormeven(n,0,1);
  x=bsxfun(@plus,x*parameters(2,:),parameters(1,:));

% uniform on [a;b]
function [x,w]=unw(n,parameters) 
  x=(1:2:2*n-1)'/(2*n);
  x=parameters(1) + (parameters(2)-parameters(1))*x;
  w=ones(n,1)/n;

% uniform on [a;b] using Gaussian quadrature nodes and weights
function [x,w]=ugnw(n,parameters) 
  [x,w]=gausslegendre(n,parameters(1),parameters(2));

  
  
function [x,w]=gnw(n,parameters)
  if size(parameters,2)>1
    error('not implemented for multiple sets of parameters')
  else
    [x,w]=qnwgamma(n,parameters(1),parameters(2));
  end
  
function [x,w]=bnw(n,parameters)
  if size(parameters,2)>1
    error('not implemented for multiple sets of parameters')
  else
    [x,w]=qnwbeta(n,parameters(1),parameters(2));
  end
  
function [x,w]=burr12nw(n,parameters)
if size(parameters,2)==1
  if n==1
    x=parameters(2)*beta(parameters(2)-1/parameters(1),1+1/parameters(1)); w=1;
  else
    % covers [0.0000001,0.9999999] of the probability range
    lb=((1-0.0000001).^(-1./parameters(2))-1).^(1./parameters(1));
    ub=((1-0.9999999).^(-1./parameters(2))-1).^(1./parameters(1));
    [x,w]=gausslegendre(n,lb,ub);
     w=w.*(parameters(1)*parameters(2)*x.^(parameters(1)-1)./(1+x.^parameters(1)).^(parameters(2)+1));
     w=w/sum(w);
  end
else
  error('not implemented for multiple sets of parameter values')
end

function [x,w]=knw(n,parameters)
if n==1
  x=parameters(2)*beta(1+1/parameters(1),parameters(2)); w=1;
else
  [x,w]=gausslegendre(n,0,1);
   a=ones(length(x),1)*parameters(1,:);
   b=ones(length(x),1)*parameters(2,:);
   xx=x*ones(1,size(a,2));
   w=bsxfun(@times,w,xx.^(a-1).*(1-xx.^a).^(b-1));
   w=bsxfun(@rdivide,w,sum(w,1));
end
  
function [x,w]=linnw(n,parameters)
if n==1
    x=0.5+parameters/6;
    w=1;
else
  [x,w]=gausslegendre(n,0,1);
  w=bsxfun(@times,w,(1-parameters)+2*x.*parameters);
end

function [x,w]=trinw(n,parameters)
if n==1
    x=sum(parameters)/3;
    w=1;
else
  a=parameters(1);
  b=parameters(2);
  c=parameters(3);
  if n>=4
    n1=max(2,min(n-2,round((c-a)/(b-a)*n)));
  elseif n==2
    n1=1;
  else
    if c>=0.5, n1=2;
    else       n1=1;
    end
  end
  [x1,w1]=gausslegendre(n1,a,c);
  [x2,w2]=gausslegendre(n-n1,c,b);
  x=[x1;x2]; w=[w1;w2];
  w=w.*((x< c).*(x-a).*(2/((b-a)*(c-a))) ...
      + (x>=c).*(b-x).*(2/((b-a)*(b-c))));
end

