% qnwgamma Quadrature nodes and weights for Gamma(lambda,theta) distribution
%   f(x) = (x/theta)^(lambda-1)*exp(-x/theta)/theta/Gamma(lambda)
%   E[x] = lambda*theta, Var[x]=lambda*theta^2 or, equivalently,
%   lambda = E[x]^2/Var[x], theta = Var[x]/E[x]
% USAGE
%   [x,w] = qnwgamma(n,lambda,theta);
%   or
%   [x,w] = qnwgamma(n,mean,variance,1);
% INPUTS
%   n      : the number of quadrature points desired
%   lambda : shape parameter of Gamma distribution [default: 1]
%   theta  : scale parmaeter of Gamma distribution [default: 1]
% The alternative syntax allows you to define the quadrature points in terms
%   of the desired mean and variance. Set the 5th input to 1 to signal that
%   the mean/variance syntax is desired.
% OUTPUTS
%   x   : prod(n) by d matrix of nodes
%   w   : prod(n) by 1 vector of weights
%
% For multivariate problems with independent Gamma random variables
% pass n and lambda as 1 by d vectors
%
% To compute expectation of f(x), where x is Gamma(lambda), write lambda
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=qnwgamma(n,lambda); E[f]=w'*f(x);
%
% Default parameter value is lambda=1 (exponential distribution)

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

function [x,w] = qnwgamma(n,lambda,theta,usemom)
d = length(n);
if nargin<2, lambda=ones(1,d); end
if nargin<3, theta=ones(1,d); end

if nargin>=4 && ~isempty(usemom) && usemom
  m=lambda; v=theta;
  lambda=m.^2./v; theta=v./m; 
end

x = cell(1,d);
[xi,w] = qnwgamma1(n(1),lambda(1));
x{1}=xi*theta(1);
for i=2:d
   [xi,wi] = qnwgamma1(n(i),lambda(i));
   x{i}=xi*theta(i);
   w=kron(w,wi); 
end
x = rectgrid(x);
return


% QNWGAMMA1 Quadrature nodes and weights for Gamma(lambda) distribution
% USAGE
%   [x,w] = qnwgamma1(n,lambda);
% INPUTS
%   n : the number of quadrature points desired
%   lambda : parameter of Gamma distribution (default=1)
% OUTPUTS
%   x : quadrature nodes
%   w : quadrature weights
% 
% To compute expectation of f(x), where x is Gamma(lambda), write lambda
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=qnwgamma(n,lambda); E[f]=w'*f(x);

function [x,w] = qnwgamma1(n,lambda)
A=diag(sqrt((1:n-1).*(lambda+(0:n-2))),-1);
A=diag(2*(1:n)+(lambda-2),0)-A-A';
[v,x]=eig(A);
[x,ind]=sort(diag(x));
w=(v(1,ind).^2)';
return

