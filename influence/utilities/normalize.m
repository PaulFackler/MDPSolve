% normalize Normalizes an array so elements of specified dimensions sum to 1
% USAGE
%    [P,ps]=normalize(P,n,sumvars);
% INPUTS
%    P         : intrinsically multidimensional array with d dimensions
%    n         : d-vector of dimension sizes (prod(n)=numel(P))
%                   [default: size(P)]
%    sumvars   : list of dimensions on which to normalize
%                   [default: 1]
%    overwrite : if P is sparse, the MEX file version of the algorithm
%                  can break MATLAB convention and overwrite P.
%                  Use with caution: non-zero value will overwrite.
%                   [default: false]
%    sp2full   : convert sparse to full if density is greater than
%                  sp2full [default: 0.5]
% OUTPUTS
%    P         : normalized array of the same size as P
%    ps        : P with dimensions in sumvars summed out
%
% This function can be used in probability calculations.
% Suppose that Pxyz=Pr(X,Y,Z|A,B,C) with elements organized
% in reverse lexicographic ordering X,Y,Z,A,B,C.
% Then [Py,Pxz]=spnormalize(Pxyz,[nx,ny,nz,na,nb,nc],2);
% will produce Py=Pr(Y|X,Z,A,B,C) and Pxz=Pr(X,Z|A,B,C);
% Reverse lexicographic ordering has first dimension changes fastest.
%
% If both X, Y and Z and A, B and C are ordered lexicographically then
%   Pxyz has dimensions [nz,ny,nx,nc,nb,na], i.e., dimensions must be
%   reversed for both rows and columns.
%
% It is assumed that ps is dense as otherwise some values of the output
% are undefined (NaN). If P is sparse these values are treated as 0 
% whereas they are set to NaN is P is full.

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

function [P,ps]=normalize(P,n,sumvars,overwrite,sp2full)
if nargin<5 || isempty(sp2full)
  sp2full=0.5;
end
if nargin<4 || isempty(overwrite)
  overwrite=false;
end
if nargin<3 || isempty(sumvars)
  sumvars=1;
end
if nargin<2 || isempty(n)
  n=size(P);
end
if prod(n)~=numel(P)
  error('prod(n) must equal numel(P)')
end
if length(size(P))>length(n)
  error('n is not compatible with the size of P')
end
if length(n)==1
  ps=sum(P);
  P=P/ps;
  return
end
if prod(n)==1
  ps=full(P);
  P=1;
  return
end
if length(n)==length(sumvars)
  if all(ismember(sumvars,1:length(n)))
    ps=sum(P(:));
    P=P/ps;
    return
  else
    error('sumvars is incorrectly specified')
  end
end
if issparse(P) && nnz(P)/numel(P)>sp2full;
  P=full(P);
end
if issparse(P) % sparse version of the algorithm
  try
    [P,ps]=spnormalizec(P,n,sumvars,overwrite);
  catch
    [P,ps]=spnormalizem(P,n,sumvars);
  end
else % full version of the algorithm
  nP=size(P);
  P=reshape(P,n);
  ps=P;
  for i=length(sumvars):-1:1
    ps=sum(ps,sumvars(i));
  end
  P=bsxfun(@rdivide,P,ps);
  P=reshape(P,nP);
  ps=squeeze(ps);
end

% m-file version of sparse algorithm
function [P,ps]=spnormalizem(P)
  [nrow,ncol]=size(P);
  cnp=[1 cumprod(n)];
  ns=n; 
  ns(sumvars)=1;
  cns=[1 cumprod(ns)];
  [ii,jj,vv]=find(P);
  indp=ii-1+(jj-1)*nrow;
  inds=1;
  for i=length(n):-1:2
    si=fix(indp./cnp(i));
    indp=indp-si*cnp(i);
    if ~ismember(i,sumvars)
      inds=inds+si*cns(i);
    end
  end
  if ~ismember(1,sumvars)
    inds=inds+indp;
  end
  ps=accumarray(inds,vv);
  vv=vv./ps(inds);
  P=sparse(ii,jj,vv,nrow,ncol);
  return
  