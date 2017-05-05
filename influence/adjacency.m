% adjacency Creates an adjacency matrix for a influence diagram
% INPUTS
%   D  : an influence diagram structure 
% OUTPUTS
%   A  : sparse adjacency matrix with A(i,j)=1 if name{i} is in the
%          parents{j} list
%   AA : sparse decendants matrix with AA(i,j)=k if there are
%          k paths from variable i to variable j

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

function [A,AA]=adjacency(D)
names=D.names;
parents=getparents(D);
n=numel(names);
if numel(parents)~=n
  error('input lists must be the same size');
end
rows=[];
cols=[];
for i=1:n
  if ~isempty(parents{i})
  [b,loc]=ismember(parents{i},1:n);
  if ~all(b)
    error(['parent list ' num2str(i) 'contains names not in list'])
  end
  rows=[rows;loc(:)]; %#ok<*AGROW>
  cols=[cols;i+zeros(numel(loc),1)];
  end
end
A=sparse(rows,cols,1,n,n);

if nargout>1
  AA=A; Ai=A; for i=1:n, Ai=Ai*A; AA=AA+Ai; end
  if any(any(Ai))
     warning('Influence diagram contains loops');
  end
end
  