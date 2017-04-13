% binprob Probability values for the binomial distribution
% Computes the probability of m successes in n trials 
%    when the success probability is p
% USAGE
%   P=binprob(p,n,m);
% INPUTS
%   p : success probability 
%   n : number of trials (non-negative integer)
%   m : number of successes (non-negative integer)
% OUTPUT
%   P : probability of m successes in n trials (same size as p)
%
% If n<m or m<0 then P=0
%
% The inputs can be scalars or conformable arrays.
% If an error occurs the procedure will attempt to break the
% problem into smaller chunks and will convert n and m to
% double.

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

function P=binprob(p,n,m)

persistent choosetab

eta=2^-1074; if eta==0, eta=realmin; end  % eta=H0000000000000001

% expand the table of choose coefficients if it is not large enough
maxnm=[max(n(:)) max(m(:))];
if any(size(choosetab)<maxnm+1)
  choosetab=gettab(maxnm(1),maxnm(2));
end

% Possible errors are out of memory or errors due to use of non-double n and/or m.
% In either case attempt to solve by looping over columns and converting to double
try
  P=(n-m).*log((1-p)+eta);
  P=exp(P+m.*log(p+eta));
  P=P.*choosetab(n+(size(choosetab,1)-1)*m+1);
catch
  if ndims(n)>2 || ndims(m)>2 || ndims(p)>2
    error('Not implemented for multidimensional arrays')
  end
  P=zeros(max(size(n),size(m)));
  if size(P,1)>size(P,2)
    nc=size(n,2);
    mc=size(m,2);
    pc=size(p,2);
    for j=1:size(P,2)
      if nc>1, nj=double(n(:,j)); else nj=double(n); end
      if mc>1, mj=double(m(:,j)); else mj=double(m); end
      if pc>1, pj=p(:,j);         else pj=p;         end
      C = nj + (size(choosetab,1)-1)*mj + 1;          % index for n choose m
      C = choosetab(min(numel(choosetab),max(1,C)));  % avoid infeasible values
      Pj=(nj-mj).*log((1-pj)+eta);
      Pj=exp(Pj + mj.*log(pj+eta));
      Pj(isnan(Pj))=1;                            % handle p(p,m) (0,0) or (1,n) 
      P(:,j)=Pj.*C;
    end
  else
    nr=size(n,1);
    mr=size(m,1);
    pr=size(p,1);
    for i=1:size(P,1)
      if nr>1, ni=double(n(i,:)); else ni=double(n); end
      if mr>1, mi=double(m(i,:)); else mi=double(m); end
      if pr>1, pi=p(i,:);         else pi=p;         end
      C = ni + (size(choosetab,1)-1)*mi + 1;          % index for n choose m
      C = choosetab(min(numel(choosetab),max(1,C)));  % avoid infeasible values
      Pi=(ni-mi).*log((1-pi)+eta);
      Pi=exp(Pi + mi.*log(pi+eta));
      Pi(isnan(Pi))=1;                            % handle p(p,m) (0,0) or (1,n) 
      P(i,:)=Pi.*C;
    end
  end
end
P((n<m) | (m<0) | (P<=exp(log(eta)))) = 0;         % set probability of infeasible outcomes to 0

% table of choose coefficients
function T=gettab(n,m)
  T=ones(n+1,m+1);
  for j=2:m+1
    T(:,j)=cumsum(T(:,j-1));
  end