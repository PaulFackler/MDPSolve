% catcountSX Creates state and state/action matrices for category count models
% USAGE
%   [S,X]=catcountSX(n,m,N,f);
% INPUTS
%   n : number of categories
%   m : number of actions
%   N : number of sites
%   f : optional handle for function indicating if an action is feasible
% OUTPUTS
%   S : n-column matrix of state variable values
%   X : n*(m+1) column matrix of state/action conbination values
%         The first n columns are the state values, the next n are the
%         values for each action associated with the sites in category 1, etc.
% S and X are returned in the smallest unsigned integer format capable of representing N.
% If these are used in arithmetic expressions they may need to be converted to
% single or double format. 
%
% f should take an n*(m+1) column matrix and return a logical vector with the same
%   number of rows as the input that indicate whether the state/action combination 
%   in a given row is feasible (1) or not (0)
% Example:
%   n=2; N=3; m=3; f=@(X) sum(X(:,n+n+1:end),2)<=1; [S,X]=catcountSX(n,m,N,f);
% Results in 
% S =
%     0     3
%     1     2
%     2     1
%     3     0
% and
% X =
%    0    3    0    2    0    0    0    1
%    0    3    0    2    0    1    0    0
%    0    3    0    3    0    0    0    0
%    1    2    0    2    0    0    1    0
%    1    2    0    2    1    0    0    0
%    1    2    1    1    0    0    0    1
%    1    2    1    1    0    1    0    0
%    1    2    1    2    0    0    0    0
%    2    1    1    1    0    0    1    0
%    2    1    1    1    1    0    0    0
%    2    1    2    0    0    0    0    1
%    2    1    2    0    0    1    0    0
%    2    1    2    1    0    0    0    0
%    3    0    2    0    0    0    1    0
%    3    0    2    0    1    0    0    0
%    3    0    3    0    0    0    0    0
% This is a situation in which there are 2 categories and 3 treatments
% but at most one site can recieve either treatment 2 or 3

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

function [S,X]=catcountSX(n,m,N,f)
% get complete set of points
if nargin<4
  S=simplexgrid(n,N,N,1);
  X=simplexgrid(n*m,N,N,1);
  X=reshape(X,[size(X,1) n m]);
  SS=zeros(size(X,1),n);
  for i=1:n
    SS(:,i)=sum(squeeze(X(:,i,:)),2);
  end
  [SS,ii]=sortrows(SS);
  X=[SS X(ii,:)];
% get only feasible points
else
  try
    [S,X]=catcountSX(n,N,m);
    feasible=logical(f(X));
    X=X(feasible,:);
  catch catcountSXError   % if out of memory error occurs build this state by state
    if ~strcmp(catcountSXError.identifier,'MATLAB:nomem')
      error(catcountSXError.message)
    end
    S=simplexgrid(n,N,N,1);
    xi=cell(1,m);
    X=[];
    ind=reshape(1:n*m,m,n)';
    for i=1:size(S,1)
      for j=1:n
        xi{j}=simplexgrid(m,double(S(i,j)),double(S(i,j)),1);
      end
      Xi=rectgrid(xi);
      Xi=[S(i+zeros(size(Xi,1),1),:) Xi(:,ind)];
      feasible=logical(f(Xi));
      X=[X;Xi(feasible,:)]; %#ok<AGROW>
    end
  end
end