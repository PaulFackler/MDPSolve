% pdfpoissonbinomial Probability Distribution Function for Poisson Binomial
% USAGE
%   P=pdfpoissonbinomial(p);
% INPUT
%   p : d-vector of success probability values
% OUTPUT
%   P : (d+1)-vector of probabilities associated with the # of successes
%
% The Poisson Binomial (or Poisson's Binomial) distribution gives the
%   probability of the number of successes in d independent trials where the 
%   probability of success on trial i is p(i). P(k) gives the
%   probability of k-1 successes.

% Uses a recursive algorithm. Let Q(i,j)=Pr(i sucesses after j trials).
% Then Q(i,j) = Q(i-1,j-1)*p(i) + Q(i,j-1)*(1-p(i))
% with appropriate endpoint adjustments.

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

function P=pdfpoissonbinomial(p)
n=length(p);
P=zeros(1,n+1);
P(1)=1-p(1);
P(2)=p(1);
for j=2:n
  pj=p(j);
  pj1=1-pj;
  P(j+1)=pj*P(j);
  for i=j:-1:2
    P(i)=pj*P(i-1)+pj1*P(i);
  end
  P(1)=pj1*P(1);
end

% P=poly((p-1)./p); P=fliplr(P/sum(P));
