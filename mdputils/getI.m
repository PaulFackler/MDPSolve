% getI Extracts information on states from state/actions combinations
% USAGE
%   [Ix,S]=getI(X,svars,S);
% or
%   [Ix,S]=getI([],svars,options);
% INPUTS
%   X     : nx x p matrix of values of the variable values
%   svars : q vector (q<=p) indicating thecolumns of X associated with the 
%             target variables
%   S     : ns x q matrix of values of the target variables
%             if omitted an index is created for the unique rows of X(:,svars)
% or
%   options : structure variable with valid fields:
%               nx     - a vector of positive integers representing the number of 
%                          values in each X variable
%                          If this is passed then X must be defined as
%                          a complete regular with prod(nx) values
%               Xint   - X is composed of integers with X(:,i) composed of values 
%                          from 0:nx(i); this allows a fast method to be used.
%               maxX   - a vector composed of the maximal values of X; used with
%                          the Xint option
%               lexico - true for lexicographic order, false for reverse lexicographic
%               scell  - second output returned as a grid cell array rather than a matrix
%               type   - determines the data type of Ix. 
%                          Valid values are 8, 16, 32 and 64 for specific unsigned integer type
%                          Also 0 for minimal unsigned integer type and -1 for double
% OUTPUTS
%   Ix    : nx-row vector of indices specifying the state associated with the
%             associated row of X
%   S     : ns x length(svars) matrix of state values
% The matrix of state values can be defined using
%   S=X(Ix,svars);
%
% Warning: if X does not contain at least one row for every value of S the
%   indices might be incorrectly computed using only 2 inputs; if this is a 
%   possability pass the complete S matrix or specify nx.
%
% getI uses one of three different methods for obtaining the index
% 1) options.nx passed and (implicitly) X is a complete regular grid with prod(nx) 
%      values such as would be created using rectgrid.
%      This is the fastest method because it avoids the construction and sorting
%      of the values matrices. Call this using
%         options=struct('nx',nx); [Ix,S]=getI([],svars,options);
%      The values of S returned are index values so S(i,j) is in {1,...,nx(1)}
% 2) Xint=true. This method assumes that X(:,1) is composed of integer values on 
%         {0,...,max(X(:,i))} (which could be on {1,...,max(X(:,i))}) 
%      It can be called using
%         options=struct('Xint',true); 
%         [Ix,S]=getI(X,svars,options);
%      or
%         options=struct('Xint',true,'maxX',max(X,[],1)); 
%         [Ix,S]=getI(X,svars,options);
%      X need not be a complete regular grid. This method requires sorting 
%         a single column of numbers.
% 3) [Ix,S]=getI(X,svars);     Uses: [S,temp,Ix]=unique(X(:,svars),'rows') 
%    [Ix,S]=getI(X,svars,S);   Uses: Ix=match(X(:,svars),S);
%    These methods are far slower than the other methods as they require
%      sorting or matching the values of each row of X(:,svars)

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011-2017, Paul L. Fackler (paul_fackler@ncsu.edu)
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

%%
function [Ix,S]=getI(X,svars,options)
if nargin<3, options=[]; S=[];
elseif isnumeric(options) && size(options,2)==length(svars)
  S=options;
  options=[];
else
  S=[];
end
nx     = [];      % sizes of X variables when X is a complete regular grid
Xint   = false;   % X(:,i) is composed of integers on {0,...,max(X(:,i))}
maxX   = [];      % max(X(:,i)); if omitted and Xint=true maxX will be calculated
lexico = true;    % lexicographic ordering is the default
type=-1;          % double is default data type
scell = false;    % if true X is returned as a grid cell array
if exist('options','var') && ~isempty(options)
  if isfield(options,'nx'),         nx         = options.nx;         end
  if isfield(options,'Xint'),       Xint       = options.Xint;       end
  if isfield(options,'maxX'),       maxX       = options.maxX;       end
  if isfield(options,'lexico'),     lexico     = options.lexico;     end
  if isfield(options,'type'),       type       = options.type;       end
  if isfield(options,'scell'),      type       = options.scell;      end
end
if ~isempty(nx),   method=1;
elseif Xint==true, method=2;
else               method=3;
end
switch method
  case 3
    if isempty(S)
      if lexico
        [S,temp,Ix]=unique(X(:,svars),'rows'); %#ok<ASGLU>
      else
        [S,temp,Ix]=unique(X(:,flipud(svars(:))),'rows'); %#ok<ASGLU>
        S=fliplr(S);
      end
    else
      Ix=match(X(:,svars),S);
    end
  case 2
    S=X(:,svars);
    if isempty(maxX)
      maxX=max(S,[],1);
    else
      maxX=maxX(svars);
    end
    [Ix,Is]=sortrowsint(S,maxX,lexico);
    if nargout>1, S=S(Is,:); end
  case 1 % nx vector passed
    if nargout>1
      [Ix,S]=getInx(nx,svars,lexico,type,scell);
    else
      Ix=getInx(nx,svars,lexico,type);
    end
    return;
end
% for methods 2 & 3 adjust output as needed
if type>=0, Ix=casttype(Ix,type); end
if nargout>1 && scell,  S=num2cell(S,1); end


%% method 1
% getI using an nx vector when X is a complete regular grid
function [Ix,S]=getInx(nx,svars,lexico,type,scell)
if iscell(nx)
  m=nx; 
  nx=cellfun(@length,m);
end
m=casttype(prod(nx(svars)),type); 

Ix=1:m;
nn=ones(1,length(nx)); nn(svars)=nx(svars);
if lexico, nn = fliplr(nn); end
if length(nn)==1; nn=[nn 1]; end  % need to pad nn due to Matlab's requirements
Ix=reshape(Ix,nn);
nn=nx; nn(svars)=1;
if lexico, nn = fliplr(nn); end
Ix=repmat(Ix,nn);
Ix=Ix(:);
if nargout>1
  S=cellfun(@(x) (1:x)',num2cell(nx(svars)),'UniformOutput',false); 
  if ~scell 
    if lexico, S=rectgrid(S);
    else       S=fliplr(rectgrid(fliplr(S)));
    end
  end
end


%% method 2
% integer valued X
function [Ix,Is]=sortrowsint(X,maxX,lexico)
if length(maxX)~=size(X,2)
  error('maxX is not compatible with X')
end
maxX=maxX(:);
if lexico
  m=flipud(cumprod([1;flipud(maxX(2:end))+1]));   % lexicographic
else
  m=cumprod([1;maxX(1:end-1)+1]);                 % reverse lexicographic
end
y=X*m;
[~,Is,Ix]=unique(y);


%% utility function
% cast a variable to a different data type
function x=casttype(x,type)
  switch type
    case  8, x=cast(x,'uint8' );
    case 16, x=cast(x,'uint16');
    case 32, x=cast(x,'uint32');
    case 64, x=cast(x,'uint64');
    case 0
      if     all(x>=2^32), x=cast(x,'uint64');
      elseif all(x>=2^16), x=cast(x,'uint32');
      elseif all(x>=2^8 ), x=cast(x,'uint16');
      else                 x=cast(x,'uint8' );
      end
    otherwise
      x=double(x);
  end