% gausslegendre Computes Guass-Legendre quadrature nodes and weights
% USAGE
%    [x,w] = gausslegendre(n,a,b);
% INPUTS
%   n   : # of nodes
%   a   : left endpoint
%   b   : right endpoint
% OUPUTS
%   x   : n by 1 vector of nodes
%   w   : n by 1 vector of weights
% 
% To compute integral of f(x) over interval [a,b], write a
% function f that returns an m-vector of values when passed an
% m by d matrix, and write [x,w]=qnwlege(n,a,b); intf=w'*f(x);
%
% Integration will be exact for polynomial of order less than 2n.
% If a and b not specified the standard interval [-1,1] is assumed.

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

function [x,w] = gausslegendre(n,a,b)
if nargin<2, a=-1; b=1; end
ii=1:(n-1);
A=diag(ii./sqrt(4*ii.^2-1),-1);
[v,x]=eig(A+A');
[x,ind]=sort(diag(x));
x=(x+1)*(b-a)/2+a;
w=(b-a)*(v(1,ind).^2)';
