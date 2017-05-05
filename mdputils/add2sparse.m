% add2sparse Places data in matrix M into a sparse matrix A
% Faciliates creating large sparse matrices by horizontal concatenation
% USAGE
%   P=add2sparse(A,M,startcol,factor);
% INPUTS 
%   A         : mxn sparse matrix
%   M         : mxq matrix
%   startcol  : starting column (scalar integer)
%   factor    : expansion factor (scalar)
%   overwrite : overwrites A in the caller's workspace
%                 is called as P=add2sparse(P,M,startcol,factor,true);
%                 P will be deleted
% OUTPUT
%   P        : mxn sparse matrix
%
% A and M must have the same number of rows
% Data from M will be placed in the q columns of A
%   starting with column startcol
% 
% Example:
%   n=50; chunk=20; nchunk=5;
%   % initialize A to be an empty sparse matrix
%   A=sparse([],[],[],n,chunck*nchunk); 
%   k=1;
%   for i=1:nchunk
%     x=sprand(n,chunk,0.5);
%     A=add2sparse(A,x,k);
%     k=k+chunk;
%   end
%
% add2sparse only uses new memory if nzmax(A) is not big enough to hold
% the new data. In this case it estimates how much more memory is required
% and creates a new sparse matrix, which is returned. The estimate of the 
% additional memory requires computes the average # of values per column
% already in the matrix (including the ones being added)times the remaining
% # of columns times an expansion factor (the default value for the
% expansion factor is 0.6).

% Uses: MEX file add2sparcec if available

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

function A=add2sparse(A,M,startcol,factor,overwrite)
if nargin<4 || isempty(factor) || factor<=0, factor=0.6; end
if nargin<5 || isempty(overwrite), overwrite=false; end
[m,n]=size(A);
[p,q]=size(M);
if (m~=p)
  error('Inputs must have the same number of rows')
end
startcol=startcol-1;
if q+startcol>n
  error('M will not fit into A')
end

% check if more memory is needed
w=nnz(A)+nnz(M);
if nzmax(A)<w
  w=ceil(w+factor*(n-startcol-q)*(w/(startcol+q))); % add estimate of needed remainder
  w=min(w,m*n);                     % can't need more than m*n
  P=sparse([],[],[],m,n,w);
  if overwrite
    A=insert(P,A,0,startcol);       % this copies old matrix to new matrix
    A=insert(A,M,startcol,q);  
  else
    P=insert(P,A,0,startcol);       % this copies old matrix to new matrix
    A=insert(P,M,startcol,q);  
  end
else
  % write M into A
  if overwrite
    A=insert(A,M,startcol,q);    % will use mex if available
  else                      
    A(:,startcol+1:startcol+q)=M;
  end
end

function A=insert(A,M,startcol,q)
try
  A=add2sparsec(A,M,startcol);  % use mex file if available
catch %#ok<CTCH>
  A(:,startcol+1:startcol+q)=M;
end