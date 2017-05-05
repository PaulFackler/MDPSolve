% sppermute Permute for an (intrinsic) sparse array 
% USAGE
%   y=sppermute(x,order,nx,nout);
% INPUTS
%   x     : sparse matrix
%   order : order for permuting the dimensions
%   nx    : intrinsic size of x
%   nout  : size of output matrix
% OUTPUT
%   y     : sparse matrix that rearranges the elements of x
%
% Equivalent to:
%   y=sparse(reshape(permute(reshape(full(x),nx),order),nout)
% but avoids the conversion from sparse to full and back
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

function y=sppermute(x,order,nx,nout)
  d=length(order);
  if any(sort(order)~=(1:d))
    error(['order must be a permutation of 1:' num2str(d)])
  end
  if nargin<3
    nx=size(x);
  else
    if numel(x)~=prod(nx)
      error('x and nx are not compatible')
    end
  end
  if length(nx)~=d 
    error('nx is not the same length as order')
  end
  if nargin<4
    nout=nx(order);
  else
    if numel(x)~=prod(nout)
      error('x and nout are not compatible')
    end
  end
  if all(diff(order)>0)   % nothing to do except reshape
    y=reshape(x,nout);
    return
  end
  if ~issparse(x)         % if x is full just use permute
    y=reshape(permute(reshape(x,nx),order),nout);
    return
  end
  
  % determine algorithm to use and call
  triedmex=false;
  try  % converting to full - generally faster for half full matrices
    if nnz(x)/numel(x)>0.5   % convert to full
      y=sparse(reshape(permute(reshape(full(x),nx),order),nout));
    else                     % try MEX function
      triedmex=true;
      y=sppermutec(x,order,nx,nout);
    end
    return
  catch
    if ~triedmex
      try % try MEX function
        y=sppermutec(x,order,nx,nout);
        return
      catch % use m-file version
        y=sppermutem(x,order,nx,nout);
      end
    else % use m-file version
      y=sppermutem(x,order,nx,nout);
    end
  end
end
    
function y=sppermutem(x,order,nx,nout)
  d=length(order);
  ny=nx(order);
  sep=find(nout(1)==[1 cumprod(ny)]);
  if isempty(sep)
    error('dimensions of y must be aligned to rows and columns of y')
  else
    sep=sep(1); % need this in case there is a singleton dimension
  end
  % get control vectors for output (y)
  yrfact=zeros(1,d);
  if sep>1
    ii=order(1:sep-1);
    yrfact(1:sep-1)=[1 cumprod(nx(ii(1:end-1)))];
  end
  ycfact=zeros(1,d);
  if sep<=d
    ii=order(sep:end);
    ycfact(sep:end)=[1 cumprod(nx(ii(1:end-1)))];
  end
  
  yrfact(order)=yrfact; ycfact(order)=ycfact;
 
  % sep is the dimension associated with the first column
  % if the input is a column vector set this to d+1
  % if the input is a row vector sep=1
  xfact=[1 cumprod(nx)];
  sep=find(size(x,1)==xfact);
  if isempty(sep)
    error('dimensions of x must be aligned to rows and columns of x')
  else
    sep=sep(1);
  end
  xfact(sep:end)=xfact(sep:end)/xfact(sep);
  
  % convert row and column values 
  [rx,cx,vx]=find(x);
  rx=rx-1;
  cx=cx-1;
  ry=ones(size(rx));
  cy=ones(size(rx));
  % work over the column dimensions of x
  if sep<=d
    for i=d:-1:sep+1
      ix=fix(cx./xfact(i));
      cx=cx-ix*xfact(i);
      if yrfact(i)>0
        ry=ry+ix*yrfact(i);
      else
        cy=cy+ix*ycfact(i);
      end
    end
    if yrfact(sep)>0
      ry=ry+cx*yrfact(sep);
    else
      cy=cy+cx*ycfact(sep);
    end
  end
  % work over the row dimensions of x
  if sep>1
    for i=sep-1:-1:2
      ix=fix(rx./xfact(i));
      rx=rx-ix*xfact(i);
      if yrfact(i)>0
        ry=ry+ix*yrfact(i);
      else
        cy=cy+ix*ycfact(i);
      end
    end
    if yrfact(1)>0
      ry=ry+rx*yrfact(1);
    else
      cy=cy+rx*ycfact(1);
    end
  end
  clear ix rx cx
  % createsparse is a faster MEX version of sparse that works when
  % all elements are unique (no adding of identically indexed values)
  y=createsparse(ry,cy,vx,nout(1),nout(2));
end

  