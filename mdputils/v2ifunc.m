% v2ifunc Return a function to find the index of the nearest neighbor on a grid
% USAGE
%   [fc,fm]=v2ifunc(x);
% INPUTS
%   x : d-element cell array with element i an n(i)x1 vector
% OUTPUT
%   fc : a function handle of the f(X1,X2,...,Xd) where Xi is an mx1 vector
%   fm : a function handle of the f(X), where X is an mxd matrix
%
% fc and fm both return an m-vector of index values on {1,..., prod(n)}
% that provides the index number associated with each element of the Xi:
%         if XX=rectgrid(x) then XX(ind,:) is the
%         nearest neighbor on the grid to X.
%
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

function [fc,fm]=v2ifunc(x)
if ~iscell(x), x={x}; end
try
  [fc,fm]=indexfunc(x);  % for evenly spaced functions this is faster
catch
  % convert x to be at the midpoint of each interval and have inf be upper bound 
  % on the last interval
  for i=1:length(x), x{i}=x{i}+[diff(x{i})/2;inf]; end
  fc=@(varargin) v2ic(x,varargin{:});
  fm=@(X)        v2im(x,X);
end

% v2i function when X is passed as separate vectors
function ind=v2ic(x,varargin)
  ind = tablookup(x{1},varargin{1});
  for i=2:d
    ind = ind*length(x{i}) + tablookup(x{i},varargin{i});
  end
  ind=ind+1;
 
% v2i function when X is passed as a single matrix
function ind=v2im(x,X)
  ind=tablookup(x{1},X(:,1));
  for i=2:d
    ind = ind*length(x{i}) + tablookup(x{i},X(:,i));
  end
  ind=ind+1;



% indexfunc Creates a function index a set of points
% USAGE
%   [Ic,Im]=indexfunc(x);
% INPUT
%   x : n x d matrix or 
%       d-element call array composed of column vectors (x <- rectgrid(x))
% OUTPUTS
%   Ic : function handle of the form I(X1,X2,...,Xd)
%   Im : function handle of the form I(X)
%
% When passed a d column vectors, each of length m, returns a vector of m integer
%   values on {1,...,n} giving the row of x associated with each of the values of X

%g={linspace(0,1,7)',linspace(1,10,10)'}; x=rectgrid(g); II=indexfunc(x); i1=randi(size(X,1),20,1); 
%X=x(i1,:); XX=mat2cell(X,size(X,1),ones(1,size(X,2))); i2=II(XX{:}); isequal(i1,i2)
function [Ic,Im]=indexfunc(x)
if iscell(x)
  x=rectgrid(x);
end
n=size(x,1);
y=(1:n)';
beta=[ones(n,1) x]\y;
alpha=beta(1);
beta=beta(2:end);
e=y-round(alpha+x*beta);
if any(e~=0)
  error('I cannot be found');
end

clear n y e x
Ic=@(varargin) Ifunc(alpha,beta,varargin{:});
Im=@(X) round(alpha+X*beta);


function ival=Ifunc(alpha,beta,varargin)
ival=alpha;
for i=1:length(beta)
  ival = ival + varargin{i}*beta(i);
end
ival=round(ival);
