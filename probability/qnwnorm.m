% qnwnorm Computes nodes and weights for multivariate normal distribution
% USAGE
%   [x,w] = qnwnorm(n,mu,var);
% INPUTS
%   n   : 1 by d vector of number of nodes for each variable
%   mu  : 1 by d mean vector
%   var : d by d positive definite covariance matrix
% OUTPUTS
%   x   : prod(n) by d matrix of evaluation nodes
%   w   : prod(n) by 1 vector of probabilities
% 
% To compute expectation of f(x), where x is N(mu,var), write a
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=qnwnorm(n,mu,var); E[f]=w'*f(x);

% Based on function of the same name in CompEcon Toolbox
% Copyright (c) 1997-2011, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu
% Modified to use eigenvalue approach

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

function [x,w] = qnwnorm(n,mu,var,tol)
d = length(n);
x = cell(1,d);
[x{1},w] = qnwnorm1(n(1));
for i=2:d
   [x{i},wi] = qnwnorm1(n(i));
   w=kron(w,wi); 
end
x = rectgrid(x);
if exist('tol','var') && ~isempty(tol)
  if tol>0
    maxw=max(w);
    ind=w>maxw*tol;
    w=w(ind); w=w/sum(w);
    x=x(ind,:);
  end
end
if exist('var','var') && ~isempty(var)
  x = x*sqrtm(var);
  if ~isreal(x)
    disp('In qnwnorm: Covariance is not positive definite - adjustments may be undesired')
    x=real(x);
  end
end
if exist('mu','var') && ~isempty(mu)
  if size(mu,1)>1, mu=mu'; end
  x=x+mu(ones(size(x,1),1),:);
end
return


% QNWNORM1 Computes nodes and weights for the univariate standard normal distribution
% USAGE
%    [x,w] = qnwnorm1(n);
% INPUTS
%   n   : number of nodes
% OUTPUTS
%   x   : n by 1 vector of evaluation nodes
%   w   : n by 1 vector of probabilities
 
function [x,w] = qnwnorm1(n)
A=diag(sqrt((1:n-1)),-1);
[v,x]=eig(A+A');
[x,ind]=sort(diag(x));
w=(v(1,ind).^2)';
return
