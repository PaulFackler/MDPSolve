% tprodm Computes products of tenors
% USAGE
%   z=tprodm(X,x2z,Y,y2z,nzout,zsparse);
% INPUTS
%   X     : a dx-dimensional tensor (implicitly)
%   x2z   : 1 or 2 row matrix with dx columns
%             x2z(1,j): the implied dimension of z associated with the jth 
%               dimension of X. Use negative values to indicate that the 
%               dimension should be summed.
%             x2z(2,j): the size (number of elements) in the jth dimension 
%               of X (if omitted it is assumed that the actual and implicit
%               sizes are equal).
%   Y     : a dy-dimensional tensor (implicitly)
%   y2z   : 1 or 2 row matrix with dy columns (analogous to x2z)
%   nzout : a vector with the desired actual sizes of z
%             (if omitted the implied sizes are used)
%   zsparse : 0/1 1 if the output is sparse [default: 1]
%               This is ignored if X and Y are both full
% OUTPUT
%   z     : the resulting tensor
%
% Tensor multiplication takes three forms:
%   1) expansion - dimension multiplies all other dimensions
%   2) matching  - dimension multiplies only matching elements in other array
%   3) summation - matches and sums over this dimension
% If both x2z and y2z contain a common element, this dimension is matched
% If the common element is negative it is also summed (which causes this
%   dimension of z to be a singleton and it can be squeezed out).
% Any value of x2z and y2z that is not in the other is an expansion dimension.
% For clarity it is best to use all unique integers for dimensions even though 
% summed diminsions are arbitrary as they are aqueezed out.
%
% For example kron(x,y) is obtained using x2z=[2 4] and y2z=[1 3].
% Ordinary matrix multiplication XY sets x2z=[1 -3] and y2z=[-3 2].
% Column-wise Kronecker products use x2z=[1 3] and y2z=[2 3].
% In each case nzout can be set if a matrix output is desired.
%
% It is not necessary for a tensor to be multi-dimensional (for sparse arrays
% only 2-D arrays are supported in MATLAB). Implicit dimensions can be used instead.
% For example, suppose that X is 20 x 15. It could implicitly be a 4-D array 
% with implicit sizes nx=[5 4 3 5]; it must be the case that prod(nx)=numel(X).

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

function z=tprodm(x,x2z,y,y2z,nzout,zsparse)
if nargin<6, zsparse=[]; end
if size(x2z,1)==2
  nx=x2z(2,:); x2z=x2z(1,:);
  %if ~issparse(x), x=reshapex(x,nx); end
else
  nx=size(x); nx(nx==1)=[];
end
if size(y2z,1)==2
  ny=y2z(2,:); y2z=y2z(1,:);
  %if ~issparse(y), y=reshapex(y,ny); end
else
  ny=size(y); ny(ny==1)=[];
end

if ~isempty(nzout) && prod(nzout)~=prod(nx(x2z>0))*prod(ny(y2z>0 & ~ismember(y2z,x2z)))
  error('nzout not compatible with other inputs')
end

% ensure that indices are unique
if length(unique(abs(x2z)))<length(x2z)
  error('x2z contains redundant indices')
end
if length(unique(abs(y2z)))<length(y2z)
  error('y2z contains redundant indices')
end

if numel(x)~=prod(nx)
  error('x and nx (x2z(2,:)) are incompatible')
end

if numel(y)~=prod(ny)
  error('y and ny (y2z(2,:)) are incompatible')
end


% combine adjacent matching subscripts
if length(x2z)>1 && length(y2z)>1 && 0
  i=1;
  while true
    if i==length(x2z), break; end
    j=find(y2z==x2z(i));
    if ~isempty(j) && j<length(y2z) && y2z(j+1)==x2z(i+1) && ...
        ~any(x2z>min(y2z(j),y2z(j+1)) & x2z<max(y2z(j),y2z(j+1))) && ...
        ~any(y2z>min(y2z(j),y2z(j+1)) & y2z<max(y2z(j),y2z(j+1))) 
      nx(i)=nx(i)*nx(i+1);
      nx(i+1)=[];
      x2z(i+1)=[];
      ny(j)=ny(j)*ny(j+1);
      ny(j+1)=[];
      y2z(j+1)=[];
    else
      i=i+1;
    end
  end
end

nz=ones(1,max(abs([x2z y2z])));
nz(abs(x2z))=nx;
nz(abs(y2z))=ny;
done=false;

if isscalar(x) 
  z=x.*y;
  [aa,bb]=sort(abs(y2z));
  z=sppermute(z,bb,ny,nzout);
  done=true;
elseif isscalar(y) 
  z=x.*y;
  [aa,bb]=sort(abs(x2z));
  z=sppermute(z,bb,nx,nzout);
  done=true;
end



% check if x has only 1 dimension
if length(x2z)==1 && any(x2z==y2z) 
  if length(y2z)==1                         % both are vectors
    if x2z<0, z=x(:)'*y(:); nz=n1;          % vectors are summed
    else      z=x.*y; nz=nx;                % vectors are matched
    end
  else
    [z,nz]=y1dim(y,y2z,ny,x,x2z,nx,nzout);  % x is vector, reverse the order
  end
  done=true;
end

% check if y has only 1 dimension
if length(y2z)==1 && any(x2z==y2z) 
  if length(x2z)==1                         % both are vectors
    if y2z<0, z=x(:)'*y(:); nz=n1;          % vectors are summed
    else      z=x.*y; nz=nx;                % vectors are matched
    end
  else
    [z,nz]=y1dim(x,x2z,nx,y,y2z,ny,nzout);  % y is vector
  end
  done=true;
end

% check if matrix methods can be used
if ~done && length(nx)==2 && length(ny)==2 
 [z,done,nz]=usemat(x,x2z,nx,y,y2z,ny,nz);
 if ~isempty(z), done=true; end
end

% if not done use fs2f, tprods or tprod
if ~done
  % check if tprod is available and use it if x and y are full
  if issparse(x) || issparse(y) || ~exist('tprod','file') 
    if ~issparse(x)  % use fs2s 
      try
        z=fs2f(x,[x2z;nx],y,[y2z;ny]);
        nz=size(z);
      catch
        [z,nz]=tprods(x,x2z,nx,y,y2z,ny,nzout,zsparse);
      end
    else
      [z,nz]=tprods(x,x2z,nx,y,y2z,ny,nzout,zsparse);
    end
  else  % use tprod MEX file
    z=tprod(reshapex(x,nx),x2z,reshapex(y,ny),y2z);
    nz=size(z); nz(nz==1)=[];
  end
end

% reshape as needed
if nargin<5 || isempty(nzout)
  if issparse(z) && length(nz)>2
     disp('cannot produce a multi-dimensional sparse array - check dimensions or set nzout')
  else
    z=reshapex(z,nz);
  end
else
  z=reshapex(z,nzout);
end

if ~isempty(zsparse)
  if zsparse
    z=sparse(z);
  else
    z=full(z);
  end
end
  
% an alternative to MATLAB's reshape that does not need the size vector
% to have at least 2 elements
function x=reshapex(x,nx)
  if length(nx)==1, x=x(:);
  else              x=reshape(x,nx);
  end
 
  
% if y has a singleton dimension pre of post multiply x with 
% either kron(I,..,y,..I) or kron(I,..,diag(y),..I) 
% (which depends on whether the matched dimension is summed or not)
function [z,nz]=y1dim(x,x2z,nx,y,y2z,ny,nzout)
  match=find(x2z==y2z);
  if isempty(match), z=[]; nz=[]; return; end
  if nx(match)~=ny
    error('matched dimensions are not of equal size')
  end
  c=cumprod(nx);
  if ~issparse(x)
    if match>1
      if c(match-1)>c(end)/c(match-1) % more rows than columns
        x=reshape(x,c(match-1),c(end)/c(match-1));
      else
        x=reshape(x,c(match),c(end)/c(match));
      end
    else
      x=reshape(x,c(1),c(end)/c(1));
    end
  end
  if size(x,1)==1 && c(1)>1
    jj=0;
  else
    jj=find(c==size(x,1));
  end

  if match<=jj  % match is in rows of x
    if match>1, m=c(match-1); else m=1; end  % size of before factors
    n=size(x,1)/c(match);                    % size of after factors  
  else          % match is in columns of x
    if match>jj+1, m=c(match-1)/size(x,1); else m=1; end   % size of before factors  
    n=c(end)/c(match);                                     % size of after factors  
  end
  if y2z>0  % no summation
    if match<=jj
      z=vxm(repmat(repmat(y(:),n,1)',m,1),x);
    else
      z=mxv(x,repmat(repmat(y(:),n,1)',m,1));
    end
  else      % sum the matching dimension
    if match<=jj
      z=kron(speye(n),kron(y(:)',speye(m)))*x;
    else
      z=x*kron(speye(n),kron(y(:),speye(m)));
    end
  end
  nz=nx(x2z>0);
  [junk,reorder]=sort(x2z(x2z>0)); %#ok<ASGLU>
  z=sppermute(z,reorder,nz,nzout);



  % see is matrix methods can be used
  % these may not be faster!
  function [z,done,nz]=usemat(x,x2z,nx,y,y2z,ny,nz)
  z=[];
  done=false;
  dz=length(nz);
  if any(size(x)~=nx)
    x=reshape(x,nx); 
  end
  if any(size(y)~=ny)
    y=reshape(y,ny); 
  end
  switch dz
    case 4 % no matches
      switch x2z(1)
        case 1
          if x2z(2)==3
            if y2z(1)==2 && y2z(2)==4
              z=kron(y,x); 
              done=true;
            elseif y2z(1)==4 && y2z(2)==2
              z=kron(y',x); 
              done=true;
            end
          end
        case 3
          if x2z(2)==1
            if y2z(1)==2 && y2z(2)==4
              z=kron(y,x'); 
              done=true;
            elseif y2z(1)==4 && y2z(2)==2
              z=kron(y',x');
              done=true;
            end
          end
        case 2
          if x2z(2)==4
            if y2z(1)==1 && y2z(2)==3
              z=kron(x,y); 
              done=true;
            elseif y2z(1)==3 && y2z(2)==1
              z=kron(x,y'); 
              done=true;
            end
          end
        case 4
          if x2z(2)==2
            if y2z(1)==1 && y2z(2)==3
              z=kron(x',y); 
              done=true;
            elseif y2z(1)==3 && y2z(2)==1
              z=kron(x',y'); 
              done=true;
            end
          end
      end
    case 3 % 1 match
      if all(x2z>0)  % no sums
        switch x2z(1)
          case 1
            switch x2z(2)
              case 2
                switch y2z(1)
                  case 1
                    z=kronrow(y,x); 
                    done=true;
                  case 3
                    z=kronrow(y',x); 
                    done=true;
                end
              case 3
                switch y2z(1)
                  case 1
                    z=kronrow(x,y); 
                    done=true;
                  case 2
                    switch y2z(2)
                      case 1
                        z=kronrow(x,y'); 
                        done=true;
                      case 3
                        z=kroncol(y,x); 
                        done=true;
                    end
                  case 3
                    z=kroncol(y',x); 
                    done=true;
                end
            end
          case 2
            switch x2z(2)
              case 1
                switch y2z(1)
                  case 1
                    z=kronrow(y,x'); 
                    done=true;
                  case 3
                    z=kronrow(y',x');
                    done=true;
                end
              case 3
                switch y2z(1)
                  case 1
                    z=kroncol(x,y);
                    done=true;
                  case 3
                    z=kroncol(x,y');
                    done=true;
                end
            end
          case 3
            switch x2z(2)
              case 1
                switch y2z(1)
                  case 1
                    z=kronrow(x',y);
                    done=true;
                  case 2
                    switch y2z(2)
                      case 1
                        z=kronrow(x',y');
                        done=true;
                      case 3
                        z=kroncol(y,x');
                        done=true;
                    end
                  case 3
                    z=kroncol(y',x');
                    done=true;
                end
              case 2
                switch y2z(1)
                  case 1
                    z=kroncol(x',y);
                    done=true;
                  case 3
                    z=kroncol(x',y');
                    done=true;
                end
            end
        end
      else           % one dimension summed (summed dimension is 3)
        switch x2z(1)
          case -3
            switch x2z(2)
              case 1
                switch y2z(1)
                  case -3
                    z=x'*y;
                    done=true;
                  case 2
                    z=x'*y';
                    done=true;
                end
              case 2
                switch y2z(1)
                  case -3
                    z=y'*x;
                    done=true;
                  case 1
                    z=y*x;
                    done=true;
                end
            end
          case 1
            switch y2z(1)
              case -3
                z=x*y;
                done=true;
              case 2
                z=x*y';
                done=true;
            end
          case 2
            switch y2z(1)
              case -3
                z=y'*x';
                done=true;
              case 1
                z=y*x'; 
                done=true;
            end
        end
      end
    case 2 % both match
      switch abs(x2z(1))
        case 1
          switch abs(y2z(1))
            case 1
              z=x.*y;
            case 2
              z=x.*y';
          end
        case 2
          switch abs(y2z(1))
            case 1
              z=x'.*y;
            case 2
              z=x'.*y';
          end
      end
      % sum as needed
      if any(x2z==-1), z=sum(z,1); end
      if any(x2z==-2), z=sum(z,2); end
      done=true;
  end
  if done
    sumdim=abs(x2z(x2z<0));
    if ~isempty(sumdim) 
      nz(sumdim)=1;
    end
  end
