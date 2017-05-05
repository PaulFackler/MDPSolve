% simplexbas Creates a basis matrix for a grid on a simplex
% A simplex is a space of q non-negative numbers that sum to C
% A regular grid on a simplex is defined by p evenly spaced intervals
%   (implying p+1 grid points for each dimension)
% USAGE
%   B=simplexbas(x,q,p,C);
% or
%   [b,ir]=simplexbas(x,q,p,C);
% INPUTS
%   x    : mx(q-1) matrix of evaluation points (each column is a point)
%            Can be mxq but the last column is ignored
%   q    : dimension of the simplex (positive integer)
%   p    : number of subintervals in each dimension (positive integer)
%   C    : grid values are non-negative and sum to C or less  (positive number)
% OUTPUTS
%   B    : nxm matrix of basis values 
%          where n is the number of grid points:
%             n = (p+q-1)!/p!/(q-1)!
% or
%   b    : qxm matrix of values
%   ir   : qxm matrix of row indices
% The two outputs are related by B=sparse(ir,ones(q,1)*(1:m),b,m,n)
%
% p=3; q=3; v=simplexgrid(q,p,p); full(simplexbas(v',q,p))
% will produce an identity matrix


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

function [b,ir]=simplexbas(x,q,p,C,checks)
if nargin<4, C=p; end
if nargin<5, checks=1; end
[m,q1] = size(x);

% error checks
if q1==q  % eliminate last column of x
  q1=q-1;
elseif q1~=q-1,
  error('x has the wrong number of columns');
end
if checks
if ~isscalar(q) || ~isscalar(C) || ~isscalar(p) 
  error('q, p and C must be scalars')
end
if (p<=0)
  error('p must be positive')
end
if (q<=1)
  error('q must be at least 2')
end
if fix(q)~=q || fix(p)~=p
  error('q and p must be integers')
end
if (C<=0)
  error('C must be positive')
end
if any(x(:)<0)
  warning('SIMPLEXBAS:Extrapolation','x contains values less than 0')
end
if any(sum(x,2)>C*(1+eps*4))
  warning('SIMPLEXBAS:Extrapolation','x contains columns that sum to more than C')
end
end

if exist('simplexbasc','file')==3  % use mex file if it exists
  if nargout<2 
    [b]=simplexbasc(x,q,p,C);
  else
    [b,ir]=simplexbasc(x,q,p,C);
  end
  return
end

% setup
p1=p+1;
if p~=C, b=(p/C)*double(x(:,1:q1)'); else b=double(x(:,1:q1)'); end
if q1==1 % faster operation when q=2 - no need to call getind because ind=v when q1=1
  v = floor(max(0,b));
  v(v>=p) = p-1; % avoids extrapolation to undefined points
  % get interpolation weights
  b = b-v;
  b = [1-b;b];  
  ir=[v+1;v+2];
else 
  tab = gettab(p1,q1);     % could pass in but it doesn't take long to compute
  % get interpolation weights
  b = flipud(cumsum(flipud(b),1));
  v = floor(b);
  v(b>=p)=p-1;
  v(b<0)=0;
  b = b-v;
  [b,ii] = sort(b,1,'descend');  % ii contains order to travel through vertices
  b = [1-b(1,:);-diff(b,1,1);b(end,:)];
  vind=(0:q1:q1*(m-1));
  ir(q,m)=0;
  ir(1,:)=getind(v,q1,p1,tab);
  for i=1:q1
    j=ii(i,:)+vind;
    v(j)=v(j)+1;
    ir(i+1,:)=getind(v,q1,p1,tab);
  end
  [ir,ii]=sort(ir,1);
  b=b(ii+ones(q,1)*(0:q:q*(m-1)));
end

% arrange output
clear v
if nargout<2
  if q1==1, n=p1; else n=tab(end,end); end
  b=sparse(ir,ones(q,1)*(1:m),b,n,m);  % convert to sparse matrix
end


function ind=getind(v,q1,p1,tab)
 ind=tab(p1+1,q1);
 q=q1+1;
 eta=p1;
 vi=v(1,:);
 for i=1:q1-1
   vi1=v(i+1,:);
   eta=eta-(vi-vi1);
   ind=ind-tab(eta,q-i);
   vi=vi1;
 end
 eta=eta-vi;
 ind=ind-tab(eta,1);
 return
 
 
% table of factorial values
function tab=gettab(p1,q1)
  tab=zeros(p1+1,q1);
  tab(:,1)=(0:p1)';
  for j=2:q1
    tab(:,j)=cumsum(tab(:,j-1));
  end
 