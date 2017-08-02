% rectbas Basis matrix for multidimensional linear interpolation
% USAGE
%   B=rectbas(x,s,evenspacing,cleanup);
% INPUTS
%   x           : mxd matrix values at which to interpolate
%   s           : a 1xd cell array of coordinate vectors
%                   each must be sorted in ascending order
%                   may be a vector rather than a cell array for d=1
%   evenspacing : 1 if all gridpoints are evenly spaced (default: 0)
%   cleanup     :  0/1/2 - Determines how extrapolation is handled
%                  0) no adjustments
%                  1) negative values are set to 0 and values are 
%                       adjusted so columns sum to 1.
%                  2) values of any variable beyond bounds are set to
%                       the nearby boundary value
% OUTPUT
%  B            : nxm interpolation matrix
%
%  n is the product of the lengths of the s vectors
%
% If the value of a function at the coordinate points is v
%   v'*B provides the linear approximation to the function at x

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

function B=rectbas(x,s,evenspacing,cleanup)
if nargin<4 || isempty(cleanup), cleanup=0; end
if isnumeric(s), s={s}; end
d=size(x,2);
if length(s)~=d
  error('s and x are incompatible in size');
end
% if necessary determine if variables are evenly spaced
if nargin<3 || isempty(evenspacing)
  evenspacing=zeros(1,d);
  for i=1:d
    ds=diff(s{i});
    if all(abs(ds-ds(1))<1e-14)
      evenspacing(i)=1;
    else
      evenspacing(i)=0;
    end
  end
elseif length(evenspacing)==1 
  evenspacing=evenspacing(ones(1,d));
elseif length(evenspacing)~=d
  error('evenspacing and s are not compatible in size')
else
  evenspacing=double(evenspacing); % rectbas1 doesn't accept anything but doubles
end

% get basis matrices
if iscell(x),  xi=x{1};
else           xi=x(:,1);
end
if cleanup==2
  xi=max(min(xi,s{1}(end)),s{1}(1));
end
B=getbas1(s{1},xi,evenspacing(1));
for i=2:d
  if iscell(x),  xi=x{i};
  else           xi=x(:,i);
  end
  if cleanup==2
    xi=max(min(xi,s{i}(end)),s{i}(1));
  end
  B=kroncol(B,getbas1(s{i},xi,evenspacing(i)));
end

if cleanup==1
  try
    B(B<0)=0;
  catch
    [ii,jj,vv]=find(B); 
    ind=vv>0;
    ii=ii(ind);
    jj=jj(ind);
    vv=vv(ind);
    B=sparse(ii,jj,vv,size(B,1),size(B,2));
  end
  B=mxv(B,1./sum(B,1),0);
  %B=bsxfun(@rdivide,B,sum(B,1));
end