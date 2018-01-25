% tablookup  Performs a table lookup
% USAGE
%   index=tablookup(tabvals,x,endadj);   returns index vector associated with x
% or
%   f=tablookup(tabvals,[],endadj);      returns function such that index=f(x)
% INPUTS
%   tabvals : a sorted vector of n values
%   x       : an array of values
%   endadj  : a optional endpoint adjustment: 0, 1, 2 or 3.
% OUTPUT
%   index   : an array (same size as x) of indices from 1 to n
%
% Returns an array of size(x) with element (i,j) equal to
%   max k: x(i,j)>=tabvals(k)
%
% Optional endpoint adjustments:
%   0: no adjustments
%   1: values of x < min(tabvals) will return 
%        length(tabvals=tabvals(1))
%   2: values of x > max(tabvals) will return 
%        m-length(tabvals=tabvals(end))
%   3: adjustments 1 and 2 will be performed
%
% With endadj=3 all the indices are between 1 and n-1
% To find the nearest table value to each x use:
%   ind = tablookup(tabvals,x,3);
%   ind = ind + (x-tabvals(ind) > tabvals(ind+1)-x);
%   nearest = tabvals(ind);

% Design and documentation based on CompEcon LOOKUP function
% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011-2017, Paul L. Fackler (paul_fackler@ncsu.edu)
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
function index=tablookup(tab,x,endadj)
if nargin<3; endadj=0; end
switch endadj
  case 0 
    tab=[tab;inf];
  case 1  
    tab(tab==tab(1))=-inf;
    tab=[tab;inf];
  case 2
    tab(tab==tab(end))=inf;
  case 3
    tab(tab==tab(1))=-inf;
    tab(tab==tab(end))=inf;
end

if isempty(x)
  clear x endadj
  index = @(x) tablookupf(x,tab);
else
  if exist('histcounts','file')  % much improved version with linear time algorithm
    [~,~,index]=histcounts(x,tab);
  else
    [~,index]=histc(x,tab);
  end
end

function index=tablookupf(x,tab)
  [~,index]=histc(x,tab);