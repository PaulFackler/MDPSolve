% rectgrid Creates a rectangular grid
% USAGE
%   X=rectgrid(x1,x2,...);
%   [X1,X2,...]=rectgrid(x1,x2,...);
%   [X{:}]=rectgrid(x1,x2,...);
% The input xi can be matrices or cell arrays of matrices 
%   (or cell arrays of cell arrays of matrices, etc.).
%   The matrices can be of any numeric data type but cannot be multi-dimensional.
% The output(s) depend on the syntax.
%   X is a matrix with every possible combination of the rows of the xi.
%     If the input matrices have different data types, the output is double.
%   [X1,X2,..] returns a set of vectors such that X=[X1,X2,...].
%   [Xc{:}] returns a cell array such that X=[Xc{:}].  Xc must be predefined as 
%     a cell array with the right number of elements. 
%   If either [X1,X2,...] or [X{:}] is used the original data types are maintained.
%     Requesting the wrong number of outputs, however, will result in an error message.
%
% Examples:
%    x1=(1:3)'; x2=(0:1)'; 
%    Xa=rectgrid(x1,x2); 
%    [Xb1,Xb2]=rectgrid(x1,x2);
%    Xc=cell(1,2); [Xc{:}]=rectgrid(x1,x2);
%         Xa =                 Xb1 = Xc{1} =        Xb2 = Xc{2} = 
%              1     0                        1                    0
%              1     1                        1                    1
%              2     0                        2                    0  
%              2     1                        2                    1
%              3     0                        3                    0 
%              3     1                        3                    1
%
%   x1=(1:3)'; x2=(0:1)'; x3=[-1;0;1]; Xa=rectgrid({x1,x2,x3}); Xb=rectgrid(rectgrid(x1,x2),x3);
%   In this example Xa and Xb are identical 18x3 matrices (so all(all(Xa==Xb)) is true).
%     
% Note: maintains lexicographic ordering.

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

function [varargout]=rectgrid(varargin)
m=1; n=0;
for i=1:nargin
  [m,n,xclass]=getcount(varargin{i},m,n,'');
end

if nargout<=1
  if strcmp(xclass,'logical')
    y=false; y(m,n)=false;
  else
    y=zeros(m,n,xclass);
  end
elseif nargout==n
  y=cell(1,n);
else
  error('Improper number of outputs requested')
end

kcol=1; k0=1; k1=m;
for i=1:nargin
  [y,kcol,k0,k1]=insert(varargin{i},y,kcol,k0,k1);
end

if nargout<=1
  varargout{1}=y;
else
  varargout=y;
end


% inserts matrices in the proper place and recursively
% handles cell arrays
function [y,kcol,k0,k1]=insert(x,y,kcol,k0,k1)
  if isempty(x), return; end
  if iscell(x)
    for i=1:numel(x)
      [y,kcol,k0,k1]=insert(x{i},y,kcol,k0,k1);
    end
  else
    m=size(x,1);
    k1=k1/m;
    ii=repmat(uint32(1:m),k1,k0);
    for i=1:size(x,2)
      if iscell(y)
        y{kcol}=x(ii,i);
      else
        y(:,kcol)=x(ii,i);
      end
      kcol=kcol+1;
    end
    k0=k0*m;
  end
    

% determines the size of the output
function [m,n,xclass]=getcount(x,m,n,xclass)
if isnumeric(x) || islogical(x)
  if ~isempty(x), m=m*size(x,1); end
  n=n+size(x,2);
  xclass=class(x);
elseif iscell(x)
  for i=1:numel(x)
    [m,n,xclass]=getcount(x{i},m,n,xclass);
  end
else
  error('Input of inproper type')
end