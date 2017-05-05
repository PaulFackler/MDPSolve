% ProjMat Forms a set of stochastic projection matrices
% USAGE
%   [A,w]=ProjMat(m,v,C,n);
% or
%   [A,w]=ProjMat(m,v,C,z,w);
% INPUTS
%   m : qx2 or qx3 matrix of mean values of (f,s) or (f,s,a)
%   v : qx2 or qx3 matrix of variance values of (f,s) or (f,s,a)
%         pass as empty or all zero matrix for a deterministic model
%   C : 2q x 2q or 3q x 3q matrix of vital rate Spearman (rank) correlations
%   n : 2q or 3q vector with number of quadrature nodes for each vital rate
%   z : N x 2q or N x 3q matrix of independent standard normal values or nodes
%   w : N-vector of associated probability weights
% OUTPUT
%   A : N-element cell array of q x q projection matrices
%   w : N-vector of associated probability weights
%
% A projection matrix is composed of fecundity rates (f), survival rates (s)
%  and, optionally, advancement rates (a)
% If there are N(i) members of class i, s(i)(1-a(i))N(i) member in stage i,
% s(i)a(i)N(i) members advance to stage i+1 and f(i)*N(i) new members of stage 1
% are produced. If m and v have 2 rows a is assumed to identically equal 1.
%
% The values of f, s, a are random and depend on correlated normal random variates z.
% f values are assumed to be Gamma variates and s and a are assumed to be Beta variates.

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

function [A,w]=ProjMat(m,v,C,n,w)
[q,d]=size(m);
if isempty(v), v=zeros(q,d); end
if ~(d==2 || d==3)
  error('m and v must have 2 or 3 columns')
end
if nargin<5  % first syntax - get quadrature nodes and weights
  if any(v(:)==0 & n(:)~=1)
    disp('In ProjMat: the number of nodes should be 1 for constant vital rates')
  end
  [z,w]=qnwnorm(n(:)');
else         % second syntax - passed quadrature nodes and weights or random values
  z=n;
end
nz=length(w);
z=z*sqrtm(2*sin((pi/6)*C));  % use transformed correlation matrix
z=real(z);
z=reshape(cdfn(z),nz,q,d);
f=zeros(nz,q);
s=zeros(nz,q);
if d==3;  a=zeros(nz,q); end
for i=1:q
  % fecundity rates
  if v(i,1)>0
    [lambda,theta]=mv2lambdatheta(m(i,1),v(i,1));
    f(:,i)=gammaincinv(z(:,i,1),lambda).*theta;
  else
    f(:,i)=m(i,1);
  end
  % survival rates
  if v(i,2)>0
    [b1,b2]=mv2ab(m(i,2),v(i,2));
    s(:,i)=icdfbeta(z(:,i,2),b1,b2);
  else
    s(:,i)=m(i,2);
  end
  % advancement rates
  if d==3 && i<q   % a(q) is identically equal to 0
    if v(i,3)>0
      [b1,b2]=mv2ab(m(i,3),v(i,3));
      a(:,i)=icdfbeta(z(:,i,3),b1,b2);
    else
      a(:,i)=m(i,3);
    end
  end
end
clear z
if d==2  % Leslie form
  A=Leslie(f,s);
else     % Lefkovitch form
  A=Lefkovitch(f,s,a);
end

% transforms mean and variance into parameters of the Gamma distribution
function [lambda,theta]=mv2lambdatheta(m,v)
lambda = m.^2./v; 
theta  = v./m;

% transforms mean and variance into parameters of the Beta distribution
function [a,b]=mv2ab(m,v)
  a=m.*(1-m)./v-1; 
  b=a.*(1-m); 
  a=a.*m; 
  
% f and s are na x q matrices of vital rates for a q-class model
% A is a cell array of na q x q projection matrices
function A=Leslie(f,s)
[na,q]=size(f);
A=cell(1,na);
for k=1:na
  A{k}=sparse([ones(1,q) 2:q q],[1:q 1:q], [f(k,:) s(k,:)],q,q);
end

% f, s and a are na x q matrices of vital rates for a q-class model
% A is a cell array of na q x q projection matrices
function A=Lefkovitch(f,s,a)
[na,q]=size(f);
A=cell(1,na);
for k=1:na
  A{k}=sparse([ones(1,q) 1:2 2:q q],[1:q 1:q 1:q], ...
              [f(k,:); s(k,:).*(1-a(k,:)); s(k,:).*a(k,:)],q,q);
end
