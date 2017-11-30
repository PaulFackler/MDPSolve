% gridmatch Finds the index of the nearest neighbor on a grid
% USAGE
%   ind=gridmatch(X,x);
% INPUTS
%   X : mxd matrix or d-element cell array composed of mx1 vectors
%   x : d-element cell array with element i an n(i)x1 vector
% OUTPUT
%   ind : mx1 vector of index values on {1,...,prod(n))
%         if XX=rectgrid(x) then XX(ind,:) is the
%         nearest neighbor on the grid to X.
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

function ind=gridmatch(X,x,evenspacing)
if ~iscell(x), x={x}; end
d=length(x);
if nargin<3 || isempty(evenspacing)
  evenspacing=false(1,d);
else
  if numel(evenspacing)==1, evenspacing=repmat(evenspacing,1,d); end
end

if iscell(X)
  ind=getind(X{1},x{1},evenspacing(1));
  for i=2:d
    ind=ind*length(x{i}) + getind(X{i},x{i},evenspacing(i));
  end
else
  ind=getind(X(:,1),x{1},evenspacing(1));
  for i=2:d
    ind=ind*length(x{i})+getind(X(:,i),x{i},evenspacing(i));
  end
end
ind=ind+1; % convert to 1-base index
    
function ind=getind(X,x,evenspacing)  
if evenspacing 
  ind=round((X-x(1))/(x(2)-x(1)));
  %if any(x(1+ind)~=X)
  %  error('oops')
  %end
else
  ind=lookup(x+[diff(x)/2;inf],X);
end
