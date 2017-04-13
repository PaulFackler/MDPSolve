% qnwnormeven Quadrature nodes and weights for Normal distribution
% Evenly spaced nodes for N(0,1) on [-a,a]
% USAGE
%   [x,w]=qnwnormeven(n,mu,sigma2);
% INPUTS
%   n      : number of nodes and weights
%   mu     : mean of the distribution
%   sigma2 : variance of the distribution
% OUTPUTS
%   x  : n-vector of nodal values
%   w  : n-vector of weights
%
% This function choses the interval that ensures
% that w'*(x-mu).^2 is as close to sigma2 as possible
% with evenly spaced nodes

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

function [x,w]=qnwnormeven(n,mu,sigma2)
if nargin<2, mu=0;    end
if nargin<3, sigma2=1; end
if n==1, x=mu; w=1; return; end
if n<=25
  options = optimset('TolFun',1e-16);
  a=fminbnd(@(a)getres(a,n),0,10,options);
  x=linspace(-a,a,n)'; 
else
  x=linspace(-8.5,8.5,n)'; 
end
w=exp(-0.5*(x.*x)); 
%w(1)=w(1)/2; w(n)=w(n)/2; % uncomment for trapezoid rule
w=w/sum(w); 
x=x*sqrt(sigma2)+mu;


function res=getres(a,n)
x=linspace(-a,a,n)'; 
w=exp(-0.5*(x.*x));
%w(1)=w(1)/2; w(n)=w(n)/2; % uncomment for trapezoid rule
w=w/sum(w); 
res=(w'*x.^2-1).^2;
