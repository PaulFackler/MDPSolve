% multisumAPPROX Normal approximation to the probabilty distribution of the sum of multinomials
% USAGE
%   [P,S]=multisumAPPROX(p,X);
% INPUTS
%   p : nxm probability matrix (p>=0 with sum(p,1)=ones(1,m))
%       Alternatively p can be a handle for a function that takes 
%         a row of X and returns an nxm probability matrix. This allows
%         p to be conditional on X.
%   X : qxm matrix of non-negative integers that sum to N
% OUTPUT
%   P : q-column matrix of probabilities with sum(P,1)=ones(1,q)
%   S : n-column matrix of values of the multinomial random variable
% The number of rows in P and S is equal to (N+n-1)!/N!/(n-1)!
%
% P(:,i) is a normal approximation to the probability distribution of S=sum(z) 
%   where z(j) is MN(p(:,j),X(i,j))
%
% NOTE: only works for q=1

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

function [P,S]=multisumAPPROX(p,X)
n=size(p,1);
p=p';
N=sum(X);
S=simplexgrid(n,N,N,1);
EX=X(:)'*p;
VX=sparse(1:n,1:n,EX,n,n)-reshape(X(:)'*kronrow(p,p),n,n);
z=double(S)-ones(size(S,1),1)*EX;
w = warning('off', 'MATLAB:nearlySingularMatrix');  %ignore nearly singular warnings
P=exp(-0.5*sum((z/VX).*z,2));
warning(w)

% set probability of impossible events to 0
elim=false;
[ii,jj]=find(p==0);
for i=1:n
  ij=ii(jj==i);
  elim=elim | (double(S(:,i))>N-sum(X(ij)));
end
[ii,jj]=find(p==1);
for i=1:n
  ij=ii(jj==i);
  elim=elim | (double(S(:,i))<sum(X(ij)));
end
P(elim)=0;
P=P/sum(P);
