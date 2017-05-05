% quasistationary Quasi-stationary distribution for an extinction process
% USAGE
%   [pi,lambda]=quasistationary(P);
% INPUTS
%   P : nxn transition matrix with one row having a single 1 on the diagonal
% OUTPUTS
%   pi     : 1x(n-1) vector of quasi-stationary distribution
%   lambda : equilibrium survival probability
%
% This function assumes that there is a single absorbing state representing
% extinction. This state is associated with the single 1 on the diagonal of P. 

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

function [pi,lambda]=quasistationary(P)
 maxit=10000;
 tol=5e-15;
 n=size(P,1);
 [ii,jj]=find(P==1);
 if length(unique(ii))~=1 
   error('not a proper extinction process')   
 end
 p=P(ii,:); p(ii)=[];
 if all(p==0)   % rows of P are current states
   pi=[0 ones(1,n-1)/(n-1)];
   for i=1:maxit
     p0=pi;
     pi=pi*P; pi(1)=0;
     pi=pi/sum(pi); 
     if max(abs(pi-p0))<tol
       break
     end
   end
   lambda=sum(pi*P);
 else       % columns of P are current states
   p=P(:,ii); p(ii)=[];
   if any(p~=0)   % rows of P
     error('not a proper extinction process') 
   end
   pi=[0;ones(n-1,1)/(n-1)];
   for i=1:maxit
     p0=pi;
     pi=P*pi; pi(ii)=0;
     pi=pi/sum(pi); 
     if max(abs(pi-p0))<tol
       break
     end
   end
   lambda=sum(P*pi);
   pi=pi';
 end
 pi(ii)=[];
