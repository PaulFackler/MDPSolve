% rvgen Generates random numbers from a specified distribution
% USAGE 
%   [x,z]=rvgen(n,rv,parents,z);
% or
%   z=rvgen(n,rv,parents,'z');
% INPUTS
%   n       : # of values to generate
%   rv      : a structure with the fields 'type' and 'parameters' (define using rvdef)
%   parents : a cell array composed of n vectors
%               If parameters is a function, parents are passed to it get
%               conditional parameter values
%   z       : n-vector of underlying random values from a previous call to rvgen
%               If the single character 'z' is passed only the underlying
%               random values are returned.
% OUTPUTS
%    x : an n-vector of random values
%    z : an n-vector of underlying random values
%          (generally these are U(0,1) except if type='n' in which case
%           they are N(0,1). These can be reused with different parameter
%           values without having to reset the generator state.
%
% See rvdef for available distributions
%
% Except type='n' values are defined using the inverse CDF method. This can
% be slow for Gamma and Beta. For Gamma consider using the Burr12, for Beta,
% consider using the Kumaraswamy.

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

function [x,z]=rvgen(n,rv,parents,z)
if nargin<3 || isempty(parents)
  parents={};
elseif isnumeric(parents)
  parents={parents};
end
% special handling for deterministic functions
if strcmp(rv.type,'f')
  % the error check really slows down the computation
  % may want to add a try/catch here to handle errors
  %if any(cellfun(@(x)any(size(x)~=[n 1]),parents)) 
  %  error('parents are not all column n-vectors')
  %end
  x=rv.valfunc(parents{:});
  z=[];
  return
end

% determine if the basic random variables (z) are needed and 
%   if that is all that is needed (getx=0)
if nargin<4, z=[]; end
getx=true;
% only get the input rvs for use in subsequent calls
if ischar(z) && strcmp(z,'z')   
  getx=false;
  z=[];
end
% get z values
if isempty(z) 
  switch rv.type
    case {'n','ne'}
      z=randn(n,1);  
    case 'i'
      z=randi(rv.parameters,n,1);
      x=z;
      return
    otherwise
      z=rand(n,1);
  end
  x=z;  % return z in the first output if getx=false
end

% get x values
if getx
  if nargin<3
    parents={};
  end
  switch rv.type
    case 'd'
      x=dgen(rv,parents,z);
    case 'logit'
      x=logitgen(n,rv,parents,z);
    case {'n','ne'}
      parameters=getparameters(n,rv,parents);
      x=z.*parameters(:,2)+parameters(:,1);
    case 'tri'
      parameters=getparameters(n,rv,parents);
      x=icdftri(z,parameters);
    case 'g'
      parameters=getparameters(n,rv,parents);
      x=gammaincinv(z,parameters(:,1)).*parameters(:,2);
    case 'b'
      parameters=getparameters(n,rv,parents);
      x=betaincinv(z,parameters(:,1),parameters(:,2));
    case 'burr12'
      parameters=getparameters(n,rv,parents);
      x=((1-z).^(-1./parameters(:,2))-1).^(1./parameters(:,1));
    case 'k'
      parameters=getparameters(n,rv,parents);
      x=(1-(1-z).^(1./parameters(:,2))).^(1./parameters(:,1));
    case 'lin'
      parameters=getparameters(n,rv,parents);
      if ~isempty(parents)
        pp=parents{1};
        a=parameters(pp,1); b=parameters(pp,2); c=parameters(pp,3); d=parameters(pp,4);
      else
        a=parameters(1,1); b=parameters(1,2); c=parameters(1,3); d=parameters(1,4);
      end
      den=a.*(d-c)+b/2.*(d.^2-c.^2);
      p1=a./den; p2=(b./den)/2; p0=-(a+b.*c/2).*c./den-z;
      x=(sqrt(p1.^2-4*p2.*p0)-p1)./p2/2;
    case 'u'
      x=z;
    case 'v'
      error('cannot simulate "v" type variables')
    otherwise
      fname=[rv.type 'gen'];
      x=feval(fname,n,rv,parents,z);
  end
  % truncate if rv is bounded
  if isfield(rv,'lb')
    x=max(x,rv.lb);
  end
  if isfield(rv,'ub')
    x=min(x,rv.ub);
  end
end

% for discrete type variables with either a single column CPT or
% with a matchfunc to match parents ro columns of the CPT
function x=dgen(rv,parents,z)
% note parameters are transposed so each row is a parameter vector
% for a different realization of the conditioning variables
if isempty(parents)
  cpt=cumsum(rv.cpt,1)';
  x = 1+sum(bsxfun(@gt,z,cpt),2); 
elseif isempty(rv.cpt)
  cpt=cumsum(rv.parameters(parents{:}))';
  x=randdisc(cpt,z);
else
  ind=rv.parameters(parents{:});
  x=randdisc(rv.cpt,z,ind);
end
if ~isempty(rv.values)
  x=rv.values(x);
end

% get the parameter values for rvs other than 'd' type
function parameters=getparameters(n,rv,parents)
if isempty(parents) || isnumeric(rv.parameters)
  parameters=rv.parameters';
else
  if ~iscell(parents) || any(cellfun(@(x)size(x,1)~=n,parents))
    if isnumeric(parents) && all(size(parents)==[n 1])
      parents={parents};
    else
      error('parents must be a cell array of n-vectors')
    end
  end
  parameters=rv.parameters(parents{:})'; % parameters are functions of parents
end

% generate values for logit random variables
function x=logitgen(n,rv,parents,z)
if isnumeric(parents)
  if isempty(parents), parents={};
  else        parents={parents};
  end
end
if ~isempty(parents)
  if size(rv.parameters,2)~=length(parents)+1
    error('parameters and # of parents are not consistent in logit type variable')
  end
  if any(cellfun(@(x)any(size(x)~=[n,1]),parents))
    error('parents must all be n-vectors');
  end
  cpt=exp([ones(n,1) [parents{:}]]*rv.parameters');
else
  cpt=exp(ones(n,1)*rv.parameters');
end
cpt=normalize(cpt);
x=randdisc(cpt,z);
if ~isempty(rv.values)
  x=rv.values(x);
end
  