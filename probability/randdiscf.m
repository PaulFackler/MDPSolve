% randdiscf Creates function for a discrete random variable simulator
% USAGE
%   f = randdiscf(p);
% INPUTS
%   p   : n x m cumulative probability matrix
% OUTPUT
%   f   : function of the form x=f(u) (if p is a vector) or x=f(u,ind) if p is a matrix)
%  
% If ind is specified
%   x(j)=i with probability p(i,ind(j))
% If ind is not specified
%   x(j)=i with probability p(i)   if p is a vector
%                           p(i,j) if q=m

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2014-2017, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function f=randdiscf(p)
if any(size(p)==1)
  c=[0;cumsum(p(:))];
  clear p
  f=@(u) randdisch1(u,c);
else
  if size(p,1)==2  % separate method to handle binary rvs
    p=p(2,:)';
    f=@(u,ind)double(u<p(ind));
  else
    if exist('randdindc','file') 
      p(end,:)=1;  % ensures that loop will stop
      f=@(u,ind)randdindc(u,p,ind);
    else
      if 1  % same algorithm as randdindc (but many times slower)
        p(end,:)=1; % ensures that loop will stop
        f=@(u,ind)randdisc(u,p,ind);
      else  % uses histc - does not appear to be as fast
        c=[zeros(1,size(p,2));cumsum(p,1)];
        c(end,:)=1;
        f=@(u,ind)randdischi(u,c,ind);
      end
    end
  end
end

 
% finds x(i) of c such that c(x(i-1),ind(i)) <= u(i) < c(x(i),ind(i)) 
% where c(i)=sum(p(1:i))
function x=randdisc(u,p,ind)
  q=length(u);
  x=zeros(size(u));
  for j=1:q
    uj=u(j);
    pj=p(:,ind(j));
    xj=0;
    while (uj>=0) 
      xj=xj+1;  
      uj=uj-pj(xj); 
    end
    x(j)=xj;
  end
  
  % finds x(i) of c such that c(x(i-1)) <= u(i) < c(x(i))
  function x=randdisch1(u,c)
    [~,x]=histc(u,c);
 
 
  % finds x(i) of c such that c(x(i-1),ind(i)) <= u(i) < c(x(i),ind(i)) 
  function x=randdischi(u,c,ind)
  q=length(u);
  x=zeros(size(u));
  for i=1:q
     [~,x(i)]=histc(u(i),c(:,ind(i)));
  end