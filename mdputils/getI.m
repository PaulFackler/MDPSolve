% getI Extracts information on states from state/actions combinations
% USAGE
%   [Ix,S]=getI(X,svars,S);
% or
%   [Ix,S]=getI(nx,svars,options);
% INPUTS
%   X     : nx-row matrix of values of the state/actions combinations; 
%              each row in a unique combination
%   svars : columns of X associated with the states
%   S     : ns x length(svars) matrix of state values
%             if omitted an index is created for the unique rows of X(:,svars)
% or
%   nx      : a vector of positive integers representing the number of values in 
%               each X variable
%   svars   : elements of X associated with the svars index
%   options : structure variable with valid fields:
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
%   possability use 3 inputs.

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

function [Ix,S]=getI(X,svars,S)
% check if X might be a positive index vactor (i.e., might actually be nx)
if isnumeric(X) && (all(size(X)>1) || any(X~=round(X)) || any(X<=0))
  if nargin<3
    [S,temp,Ix]=unique(X(:,svars),'rows'); %#ok<ASGLU>
  else
    Ix=match(X(:,svars),S);
  end
else % nx vector passed
  nx=X; 
  if nargin>2, options=S; end
  lexico = true;    % lexicographic ordering is the default
  type=-1;          % double is default data type
  scell = false;    % if true X is returned as a grid cell array
  if exist('options','var') && ~isempty(options)
    if isfield(options,'lexico'),     lexico     = options.lexico;     end
    if isfield(options,'type'),       type       = options.type;       end
    if isfield(options,'scell'),      type       = options.scell;      end
  end
  if iscell(nx)
    m=nx; 
    nx=cellfun(@length,m);
  end
  m=prod(nx(svars));

  % determine the data type of the index
  switch type
    case 8,  m=cast(m,'uint8');
    case 16, m=cast(m,'uint16');
    case 32, m=cast(m,'uint32');
    case 64, m=cast(m,'uint64');
    case 0
      if     m>=2^32, m=cast(m,'uint64');
      elseif m>=2^16, m=cast(m,'uint32');
      elseif m>=2^8,  m=cast(m,'uint16');
      else            m=cast(m,'uint8');
      end
    otherwise
      m=double(m);
  end
  Ix=1:m;

  nn=ones(1,length(nx)); nn(svars)=nx(svars);
  if lexico, nn = fliplr(nn); end
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
end