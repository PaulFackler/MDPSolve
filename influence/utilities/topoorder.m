% topoorder Topological ordering for a Directed Acyclic Graph (DAG)
% USAGE
%   [L,errflag]=topoorder(A);
% INPUT
%   A : an adjacency matrix or an influence diagram structure
% OUTPUTS
%   L : a vector with a variable ordering
%   errflag : 0-successful 1-graph not a DAG (returns empty L)

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

function [L,errflag]=topoorder(A)
if isstruct(A)
  A=adjacency(A);
end
errflag=0;
n=size(A,1);
k=0;
L=zeros(1,n);
S=find(sum(A,1)==0);
while 1
  if isempty(S), break; end
  next=S(1); S(1)=[];
  k=k+1; L(k)=next;
  children=find(A(next,:)~=0);
  A(next,:)=0;
  for i=1:length(children)
    if all(A(:,children(i))==0)
      S=[S children(i)];
    end
  end
end
if any(any(A))
  L=[];
  if nargout<2
    warning('graph is not acyclic'); 
  else
    errflag=1;
  end
end
