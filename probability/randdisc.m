% randdisc Discrete random variable simulator
% USAGE
%   x = randdisc(p,u,ind);
% INPUTS
%   p   : n x m probability matrix
%   u   : q x 1 vector of random uniform variates
%   ind : q x 1 vector of indices for columns of p (not used if m=1)
% OUTPUT
%   x   : q x 1 vector of random values on {1,...,n}
%  
% If ind is specified
%   x(j)=i with probability p(i,ind(j))
% If ind is not specified
%   x(j)=i with probability p(i)   if p is a vector
%                           p(i,j) if q=m

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

function x=randdisc(p,u,ind)
if any(size(p)==1)
  cp=[0;cumsum(p(:))];
  [~,x]=histc(u,cp);
  return
else
  q=numel(u);
  if nargin<3
    if size(p,2)~=q
      error('cannot tell which columns of p to use')
    end
    ind=[];
  else
    if numel(ind)~=q
      error('ind must be the same size and u')
    end
  end
  if size(p,1)==2  % separate method to handle binary rvs
    if nargin<3
      x=1+double(u>p(1,:)');
    else
      x=1+double(u>p(1,ind)');
    end
    return
  end
  %p=cumsum(p,1); p=[zeros(1,size(p,2));p(1:end-1,:)]; 
  %if nargin<3
  %  x=sum(bsxfun(@ge,u(:)',p),1)';
  %else
  %  x=sum(bsxfun(@ge,u(:)',p(:,ind)),1)';
  %end
  %return
  %try
  %  error(' ')
  %  if isempty(ind)
  %    x=randdiscc(p,u,1:q);
  %  else
  %    x=randdiscc(p,u,ind);
  %  end
  %catch
    m=size(p,1);
    x=zeros(size(u));
    for j=1:q
      z=u(j);
      if nargin<3
        pj=p(:,j);
      else
        pj=p(:,ind(j));
      end
      for i=1:m;  
        z=z-pj(i);  
        if z<0
          break;
        end
      end
      x(j)=i;
    end
  %end
end