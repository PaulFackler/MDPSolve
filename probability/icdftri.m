% icdftri Inverse CDF of the triangular distribution
% USAGE
%   x=icdftri(u,a,b,c);
% or
%   x=icdftri(u,parameters);
% INPUTS
%   u : vector of values on [0,1]
%   a : lower bound
%   b : uppser bound
%   c : model value 
% OUTPUT
%   x : vector of values on [a,b]
%
% Note: it must be the case that a<=c<=b (not checks are made)
%
% a, b and c can be scalars or vectors of the same size as u
%
% Can also be called as 
%    x=icdftri(u,[a,b,c]);
% or 
%    x=icdftri(u,{a,b,c});

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

function x=icdftri(u,a,b,c)
if nargin==2
  if iscell(a)
    c=a{3};
    b=a{2};
    a=a{1};
  else
    c=a(:,3);
    b=a(:,2);
    a=a(:,1);
  end
end
ba=b-a;
bc=b-c;
ca=c-a;
ind=u<ca./ba;
x=(a+sqrt(u.*ba.*ca)).*ind + (b-sqrt((1-u).*ba.*bc)).*(~ind);