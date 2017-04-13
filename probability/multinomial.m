% multinomial Probability of a multinomial random variable
% USAGE
%   prob=multinomial(p,n)
% INPUTS
%   p    : m-vector of probabilities (non-negative values that sum to 1)
%   n    : m-vector of results of N=sum(n) trials
% OUTPUT
%   prob : probability of obtaining n

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

function prob=multinomial(p,n)
if abs(sum(p)-1)>1e-15
  warning('Multinomial:sum1','p may not sum to 1')
end
m=length(p);
if numel(n)==1
  N=n;
  n=simplexgrid(m,N,N);
  prob=N;
else
  if size(n,2)~=m
    if size(n,1)==m && size(n,2)==1
      n=n';
    else
      error('n is incompatible with p')
    end
  end
  prob=sum(n,2);
  N=max(prob);
end
gln=log([1 cumprod(1:N)]);
eta=2^-1074; if eta==0, eta=realmin; end  % eta=H0000000000000001
prob=exp(n*log(p(:)+eta) + gln(prob+1)' - sum(gln(n+1),2));
prob(prob<=exp(log(eta)))=0;
