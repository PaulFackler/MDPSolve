% simplexindex Converts simplex vertices into index values
% USAGE
%   index=simplexindex(v,q,p,C,tab);
% or
%   tab=simplexindex([],q,p,C);
% INPUTS
%   v     : mx(q-1) matrix of evaluation points
%             Can be mxq but the last column is ignored.
%   q     : dimension of the simplex (positive integer)
%   p     : number of subintervals in each dimension (positive integer)
%   C     : grid values are non-negative and sum to C or less (positive number)
%             [default: p]
%   tab   : table of multiset coefficients. For repeated calls to simplexindex
%             obtain tab using the second syntax and pass it on subsequent calls.
% OUTPUT
%   index : mx1 vector of index values

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

function index=simplexindex(v,q,p,C,tab)
if nargin<3, error('3 inputs are required'); end
p1=p+1;
q1=q-1;
if nargin<4, C=p; end
if nargin<5 || isempty(tab)
  tab=gettab(p1,q1);
  if isempty(v)
    index=tab;
    return
  end
end
%try
%  index=simplexindexc(v,q,p,C,tab);
%catch
  if C==p
    index=getindI(v,q1,p1,tab);
  else
    index=getindD(v,q1,p1,tab,double(p)/double(C));
  end
%end
end


function ind=getindI(v,q1,p1,tab)
  ind=tab(p1+1,q1);
  datatype=class(ind);
  q=q1+1;
  eta=p1;
  for i=1:q1
    eta=eta-cast(v(:,i),datatype);
    ind=ind-tab(eta,q-i);
  end
end
 
 
 function ind=getindD(v,q1,p1,tab,factor)
   ind=tab(p1+1,q1);
   q=q1+1;
   eta=p1;
   for i=1:q1
     eta=eta-round(v(:,i)*factor);
     ind=ind-tab(eta,q-i);
   end
 end

 
% table of factorial values
function tab=gettab(p1,q1)
  tab=[zeros(1,q1);ones(p1,q1)];
  tab(:,1)=(0:p1)';
  for j=2:q1
    tab(:,j)=cumsum(tab(:,j-1));
  end
  return
  [ct,maxsize]=computer;
  if maxsize==2147483647  % 32 bit machine
    tab=uint32(tab);
  else                    % 64 bit machine
    tab=uint64(tab);
  end
end