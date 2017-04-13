% mfpt Mean First Passage Time for Markov chains
% USAGE
%   [M,w]=mfpt(P);
% or
%   M=mfpt(P,target);
% INPUTS
%   P      : nxn transition probability matrix (non-negative element w/ row sums = 1)
%   target : an n-element logical vector or a vector of unique integer values in {1,...,n}
%              This vector specifies the target set of states.
% OUTPUTS
%   M : mean first passage time matrix: M(i,j) is the mean time to first reach 
%          state j starting at state i
%          or
%          n-vector where M(i) is the mean time to first reach one of the target 
%          states starting at state i
%   w : long-run distribution
%
% Note: with the first syntax a check is made for recurrence and a warning is issued
%   If the chain is not recurrent, the method used to determine M may produce
%   incorrect results.
%
% P can be a column stochastic matrix (columns sum to 1) but if it is doubly 
%   stochastic (rows and columns both sum to 1) mfpt interprets rows as being 
%   associated with current states but prints a warning message. If P is doubly
%   stochastic with columns associated with current states pass in P'.

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

function [M,w]=mfpt(P,target)
  n=size(P,1);
  if any(abs(sum(P,2)-1)>1e-15*n)
    if any(abs(sum(P,1)-1)>1e-15*n)
      error('probability matrix appears to be invalid')
    end
    P=P';
  else
    if all(abs(sum(P,1)-1)<1e-15*n)
      disp('In mfpt: P appears to be doubly stochastic - assuming rows are current')
    end
  end
  if nargin<2
	M = inv(full(P)-eye(n)+1);
	w = sum(M,1)';
  if any(w<1e-8/n)
    disp('Probability matrix is not recurrent')
  end
  for c = 1:n
		M(:,c) = (M(:,c) - M(c,c))./w(c);
  end
  else
    if islogical(target)
      if length(target)~=n
        error('if target is logical it must contain n elements')
      end
      n0=n-sum(target);
    else
      n0=n-length(target);
      if any(target<1 || target>n)
        error('target vector is incorrectly specified');
      end
      ind=target;
      target= false(n,1);
      target(ind)=true;
    end
    if n0<=0
      error('target vector is incorrectly specified');
    end
    p=sum(P(~target,target),2);
    if all(p==0)   % target set is inaccessible from initial set
      M=zeros(n,1)+inf;
    else
      PP=[2 sparse(1,n0);
          p P(~target,~target)];
      M=zeros(n,1);
      temp=(speye(1+n0)-PP)\[0;ones(n0,1)];
      M(~target)=temp(2:end);
    end
  end