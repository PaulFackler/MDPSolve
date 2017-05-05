% f2P Creates a discrete transition probability matrix from a conditional density
% USAGE
%   P=f2P(f,S,X,w,options);
% INPUTS
%   f : density function of form f(S,X,...)
%   S : ns-row matrix of state variable values
%   X : nx-row matrix of state/action combinations
%   w : ns-vector of quadrature weights [default: w=1]
%   options: a structure variable
%           tol: tolerance for setting probabilities to 0
%       can also pass a value for tol (0<=tol<small value)
% OUTPUTS
%   P : ns x nx matrix of transition probabilities (column stochastic)
%
% f is the conditional density of S+ given X. Its syntax is
%   f(Splus,X)
% It should accept an nsxd matrix S and a 1xq vector X and return an
%   ns-vector of non-negative probability density values.
%
% If tol>0 the transition matrix is sparse with values less than tol/ns
% set to 0.
%
% This procedure uses the idea that
%   E[V+|X] = int f(S+|X)V(S+)dS+ 
% is approximately equal to 
%    sum w(i)*f(S(i)|X))*V(S(i))
% Thus the probabilities Prob(S+|X) can be approximated using
%   P(S(i)|X(j)) = w(i)f(S(i)|X(j))

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function P=f2P(f,S,X,w,options)
% get 0-tolerence level
if iscell(S)
    S=rectgrid(S);
end
tol = 0; 
if exist('options','var') && ~isempty(options)
  if isnumeric(options)
    tol=options/size(S,1);
  else
    if isfield(options,'tol');
      tol=options.tol/size(S,1);
    end
    if isfield(options,'byS');
      byS=options.byS;
    end
  end
end
ns=size(S,1);
nx=size(X,1);
% check if quadrature weights were passed
if exist('w','var') && ~isempty(w)
  hasweights=true;
  if numel(w)~=ns
    error('Number of rows of S must equal the length of weights vector w')
  end
  w=w(:);
else
  hasweights=false;
end
S=S+realmin;
X=X+realmin;
% if tol>0 get a sparse transition matrix
if tol>0
  P=sparse(ns,nx);
  for j=1:nx
    Pj=f(S,X(j,:));
    Pj(Pj<tol*max(Pj(:)))=0;
    P=add2sparse(P,sparse(Pj(:)),j,0.6,true);
  end
% if tol=0 get a dense transition matrix
else
  P=zeros(ns,nx);
  byS=1;
  if byS
    for i=1:ns
      P(i,:)=f(S(i,:),X);
    end
  else
    for j=1:nx
      P(:,j)=f(S,X(j,:));
    end
  end
end
if hasweights, P=vxm(w,P); end
% normalize to sum to 1
Psum=sum(P,1);
if any(isnan(Psum))
  disp(['In f2P:  number of columns containing invalid values: '  num2str(sum(isnan(Psum)))])
  Psum(isnan(Psum))=1;
end
if any(Psum==0)
  disp(['In f2P:  number of columns summing to 0: '  num2str(sum(Psum==0))])
  Psum(Psum==0)=1;
end
P=mxv(P,1./Psum);
