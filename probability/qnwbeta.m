% qnwbeta Computes quadrature nodes and weights for Beta(a,b) distribution
% USAGE
%   [x,w] = qnwbeta(n,a,b);
% INPUTS
%   n   : number of quadrature points desired
%   a,b : Beta distribution parameters  
% OUTPUTS
%   x   : prod(n) by d matrix of nodes
%   w   : prod(n) by 1 vector of weights
%
% For multivariate problems with independent Beta random variables
% pass n, a, and b as 1 by d vectors
%
% To compute expectation of f(x), where x is Beta(a,b), write a
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=qnwbeta(n,a,b); E[f]=w'*f(x);
%
% Default parameter values are a=b=1 (uniform distribution)

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
 
function [x,w] = qnwbeta(n,a,b)

d = length(n);
if nargin<2, a=ones(1,d); end
if nargin<3, b=ones(1,d); end

x = cell(1,d);
[x{1},w] = qnwbeta1(n(1),a(1),b(1));
for i=2:d
   [x{i},wi] = qnwbeta1(n(i),a(i),b(i));
   w=kron(w,wi); 
end
x = rectgrid(x);
return


% QNWBETA1  Computes quadrature nodes and weights for Beta(a,b) distribution
% USAGE
%   [x,w] = qnwbeta1(n,a,b);
% INPUTS
%   n   : number of quadrature points desired
%   a,b : Beta distribution parameters
% OUTPUTS
%   x   : quadrature nodes
%   w   : quadrature weights
function [x,w]=qnwbeta1(n,a,b)
a=a-1;
b=b-1;
i=1:n;
ii=2*i+a+b;
s = warning('off','MATLAB:divideByZero');
if a==b
  aa=ones(n,1)/2;
else
  aa=(1-(b^2-a^2)./(ii-2)./ii)/2;
end

i(n)=[];
ii(n)=[];
bb=sqrt(i.*(i+a).*(i+b)./(ii+1)).*sqrt((i+a+b)./(ii-1))./ii;
warning(s)  % restore warning state

if abs(a+b+1)<2*eps,  aa(1)=(1+a-b)/2; bb(1)=sqrt(a*b/2);          end
if abs(a+b)  <2*eps,  aa(1)=(2+a-b)/4; bb(1)=sqrt((a+1)*(b+1)/12); end

A=diag(bb,-1);
A=diag(aa)-A-A';
[v,x]=eig(A);
[x,ind]=sort(diag(x));
w=(v(1,ind).^2)';

return


