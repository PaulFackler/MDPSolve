% freudenthal Freudenthal interpolation basis matrix for a lattice
% USAGE
%   B=freudenthal(S,s,cleanup);
% INPUTS
%   S    : Nxd matrix of evaluation points
%   s    : d element cell array of vectors defining the lattice
%   cleanup :  0/1/2 - Determines how extrapolation is handled
%              0) no adjustments
%              1) negative values are set to 0 and values are 
%                   adjusted so columns sum to 1.
%              2) values of any variable beyond bounds are set to
%                   the nearby boundary value
% OUTPUT
%   B    : nxN matrix of basis values (n is the number of grid points)
%
% s defines a regular lattice containing n points, where n is the
% product of the numbers of values in vectors in s. If f is an 1xn vector
% of function values at the n points, then f*B is a 1xN vector of 
% interpolated values of the function, evaluated at X. 
%
% Values of B will always be non-negative if all evaluation points are
% inside the grid. Columns of B sum to 1.
%
% This implements what is sometimes called Coxeter-Freudenthal-Kuhn 
% triangulation.
%
% Note that the rows of B are arranged so B'*rectgrid(s) equals S.

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

function B=freudenthal(S,s,cleanup)
if nargin<3, cleanup=false; end
if iscell(s)
  d=length(s);
else
  d=1;
  s={s};
end
[N,dim]=size(S);
if dim~=d
  error('s and S are incompatible')
end
  
n=zeros(dim,1);
ind=zeros(N,dim);
d=zeros(N,dim);
for i=1:dim
  n(i)=length(s{i});
  df=diff(s{i});
  Si=S(:,i);
  if cleanup==2 
    Si=max(min(Si,s{i}(end)),s{i}(1));
  end
  if all(abs(df./df(1)-1)<1e-15)
    hi=df(1);
    ind(:,i)=evenlookup(s{i},Si,hi);
  else
    ind(:,i)=lookup(s{i},Si,3);
    hi=s{i}(ind(:,i)+1)-s{i}(ind(:,i));
  end
  d(:,i)=(Si-s{i}(ind(:,i)))./hi;
end

[d,p]=sort(d,2,'descend');

% compute weights
lambda=[-diff(d,1,2) d(:,end)];
lambda=[1-sum(lambda,2) lambda];
if cleanup==1
  lambda(lambda<0)=0;
  lambda=vxm(1./sum(lambda,2),lambda);
end

% get the row indices for B
bw=[flipud(cumprod(flipud(n(2:end))));1];
ind=(ind-1)*bw+1;
bind=zeros(N,dim+1);
bind(:,1)=ind;
for i=1:dim
  ind=ind+bw(p(:,i));
  bind(:,i+1)=ind;
end
B=sparse(bind,(1:N)'*ones(1,dim+1),lambda,prod(n),N);
% uncomment next line (and comment previous one) to transpose B
%B=sparse((1:N)'*ones(1,dim+1),bind,lambda,N,prod(n));
end

% Finds the low vertex of the simplex containing S
% when s is evenly spaced
function ind=evenlookup(s,S,h)
  n=length(s);
  ind=ceil((S-s(1))/h);
  ind=min(max(ind,1),n-1);
end

