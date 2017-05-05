% match Finds the index value of the closest element in a point set
% USAGE
%   ind = match(S,s);
% INPUTS
%   S   : a pxd matrix
%   s   : an nxd matrix
% OUTPUT
%   ind : a p-vector of integers in {1,...,n} indicating the row of
%         s that most closely matches each row in S (using Euclidean norm)
%
% Example:
%  s=[0 0; 0 1; 1 0; 1 1]; S=[0 0; 0 0; 1 0; 0 0; 1 1];
%  match(s,S) returns
%  [1; 1; 3; 1; 4]

% The kdtree and kdtree_closestpoint mex code is due to Steven Michael
% http://www.mathworks.com/matlabcentral/fileexchange/7030-kd-tree-nearest-neighbor-and-range-search
% Copyright (c) 2005, Steven Michael
% All rights reserved.

% The matchm code is an edited version of the function knnsearch by Yi Cao
% By Yi Cao at Cranfield University on 25 March 2008
% Copyright (c) 2009, Yi Cao
% All rights reserved.
% http://www.mathworks.com/matlabcentral/fileexchange/19345-efficient-k-nearest-neighbor-search-using-jit

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

function ind = match(S,s)
if isfloat(S) && isfloat(s)
  try
    tree = kdtree(s);
    ind  = kdtree_closestpoint(tree, S);
  catch
    ind = matchm(S,s);
  end
else
  ind = matchm(S,s);
end


function [ind,dist]=matchm(S,s)
[m,n] = size(S);
p=size(s,1);
ind = zeros(m,1);
dist = 0;
for i=1:m
  d=zeros(p,1);
  for j=1:n
    d=d+(double(s(:,j))-double(S(i,j))).^2;
  end
  [dist,ind(i)]=min(d);
end

