% simplexgrid Creates a grid for a simplex of dimension q
% Uses p subintervals spanning the integers 0 to p
% USAGE
%   S=simplexgrid(q,p,C,getall,type,Smax);
% INPUTS
%   q       : scalar problem dimension or number of variables (positive integer)
%   p       : scalar number of subdivisions for each dimension (positive integer)
%   C       : grid points sum to C or less (default: p)
%   getall  : returns an nxq matrix with each row summing to exactly C [default]
%             if getall==0 only q-1 columns are returned
%   type    : string specifying the data type to use (uint8, uint16, uint32, double)
%             alternatively: 1 for double [the default], 
%                            0 to use the smallest possible unsigned integer
%             note: if C~=p x is double and type is ignored
%   Smax    : q-vector of maximal values of each column of S
% OUTPUTS
%   S  : nxq matrix of values with each row summing to C 
%            where n=(p+q-1)!/p!/(q-1)! (the multiset(p,q) coefficient)
%            If Smax is set n can be far smaller.
%        (alternate: use getall=false to produce nx(q-1) matrix with rows
%            summing to C or less)
%
% If C=p the grid values can be returned as unsigned integers. This saves space
% but limits the usefulness of the grid for doing arithmetic. 

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

function x=simplexgrid(q,p,C,getall,type,smax)

if ~exist('C','var')       || isempty(C),       C=p;          end
if ~exist('getall','var')  || isempty(getall),  getall=true;  end
if ~exist('type','var')    || isempty(type),    type=1;       end
if ~exist('Smax','var')    || isempty(smax),    smax=[];      end

if isnumeric(type) && type~=0
    type='double';
end
if isnumeric(type)
  if C~=p
    type='double';
  elseif p<256
    type='uint8';
  elseif p<65536
    type='uint16';
  elseif p<4294967296
    type='uint32';
  else
    error('problem too big');
  end
end
if nargin<6 || isempty(smax) || all(smax>=p)
  x=sgnomax(q,p,type);
else
  if length(smax)~=q
    error('smax must be a q-vector')
  end
  if sum(smax)<p
    x=zeros(0,q);  % no feasible values
    return
  end
  x=sgwithmax(q,p,type,smax);
end
if ~getall
  x=x(:,1:q-1);
end
if C~=p
  x=x*(double(C)/double(p));
end
end


function x=sgnomax(q,p,type)
m=multiset(q,p);
x=zeros(m,q,type);
q1=q-1;
xn=[zeros(1,q1,type) p];
n=1;
x(1,:)=xn;
while n<m
  if xn(q)>0
    xn(q1)=xn(q1)+1;
    xn(q)=xn(q)-1; 
  else
    k=q1;
    while xn(k)==0, k=k-1; end
    xn(k-1)=xn(k-1)+1;
    xn(q)=xn(k)-1;
    xn(k)=0;
  end   
  n=n+1;
  x(n,:)=xn;
end
end
  

function x=sgwithmax(q,p,type,smax)
m=multiset(q,p);
try
  x=zeros(m,q,type);
catch SIMPLEXGRIDerror
  if ~strcmp(SIMPLEXGRIDerror.identifier,'MATLAB:nomem')
    rethrow(SIMPLEXGRIDerror)
  end
  m=sgwithmaxNum(q,p,type,smax);  % determine size of output
  x=zeros(m,q,type);
end  
q1=q-1;
xn=[zeros(1,q1,type) p];
k=q1;
n=0;
smax1=smax(1);
smaxq=smax(q);
while (xn(1)<=smax1)
  while true
    while xn(k)>smax(k)
      if k==1; break; end
      xn(k-1)=xn(k-1)+1;
      xn(q)=xn(q)+xn(k)-1;
      xn(k)=0;
      k=k-1;
    end
    if xn(q)>smaxq
      xn(q1)=xn(q)-smaxq;
      xn(q)=smaxq;
      k=q1;
    else
      break;
    end
  end
  if xn(1)>smax1, break; end
  n=n+1;
  x(n,:)=xn;
  k=q1;
  if xn(q)>0
    xn(k)=xn(k)+1;
    xn(q)=xn(q)-1;
  else
    while xn(k)==0, k=k-1; end
    if k==1, break; end
    xn(k-1)=xn(k-1)+1;
    xn(q)=xn(k)-1;
    xn(k)=0; 
    k=k-1;
  end
end
x=x(1:n,:);
end


function n=sgwithmaxNum(q,p,type,smax)
q1=q-1;
smaxq=smax(q);
smax1=smax(1);
xn=zeros(1,q1,type);
xnq=p;
k=q1;
n=0;
while (xn(1)<=smax1)
  while true
    while xn(k)>smax(k)
      if k==1; break; end
      xn(k-1)=xn(k-1)+1;
      xnq=xnq+xn(k)-1;
      xn(k)=0;
      k=k-1;
    end
    if xnq>smaxq
      xn(q1)=xnq-smaxq;
      xnq=smaxq;
      k=q1;
    else
      break;
    end
  end
  if xn(1)>smax1, break; end
  n=n+1;
  k=q1;
  if xnq>0
    xn(k)=xn(k)+1;
    xnq=xnq-1;
  else
    while xn(k)==0, k=k-1; end
    xn(k-1)=xn(k-1)+1;
    xnq=xn(k)-1;
    xn(k)=0; 
    k=k-1;
  end
end
end

function m=multiset(q,p)
  m=1:p+1;
  for i=1:q-2
    m=cumsum(m);
  end
  m=m(end);
end
