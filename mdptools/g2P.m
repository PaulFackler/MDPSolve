% g2P Creates a discrete transition probability matrix
% USAGE
%   [P,warnings]=g2P(g,s,X,e,w,options);
% INPUTS
%   g : function of form g(X) or g(X,e)
%   s : grid information about state variables 
%         (a cell array of state variable vectors)
%         The matrix of state values (S) should be obtainable
%         using S=rectgrid(s);
%   X : k-row matrix of state/action combinations
%   e : pxq matrix of random shocks [optional]
%   w : p-vector or kxp matrix of probability weights [optional]
%         Alternative: pass e and w as m-element cell arrays of vectors
%           which will be combined using disccombine (this assumes that
%           the m-noise variables are independent)
%   options : a structure variable to control the procedure
%               cleanup    : 0/1/2 - Determines how extrapolation is handled
%                            0) no adjustments
%                            1) negative values are set to 0 and values are 
%                                 adjusted so columns sum to 1.
%                            2) values of any variable beyond bounds are set to
%                                 the nearby boundary value
%               rectinterp : 0) use Freudenthal interpolation 
%                            1) use linear interpolation [default: 1]
%               expande    : 0) e is passed to g as a row vector
%                            1) e is passed as a kxq matrix with equal rows
%                            [default: 0]
%             Alternatively a scalar (0/1/2) can be passed for the value of cleanup
%             and rectinterp defaults to 1.
% OUTPUTS
%   P        : nxk matrix of transition probabilities (column stochastic)
%   warnings : a cell array of warning messages (if not requested any messages
%                will be displayed)
%
% The function g should accept a k-row matrix of state/action combinations
%   (each row with a different combination) and return a k-row matrix
%   of next period state values. If the problem is stochastic (e is
%   non-empty) g should accept a k-row matrix X and an associated k-row
%   matrix of random noise terms.
%
% Currently this function only handles regular rectangular grids.

% programming note:
% the function implements interpolation for simplex and scattered grids as well 
% but these features have not be adequately tested or documented

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011-2013, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function [P,warnings]=g2P(g,s,X,e,w,options)
warnings={};
if nargin<6, options=0; end
geometry=0;
if isnumeric(options)
  cleanup=options;
  rectinterp=1;
  expande=false;
else
  getopts(options, ...
         'cleanup',     0, ... % ensures that P is a probability matrix
         'rectinterp',  1, ... 
         'expande', false  ); 
end
if rectinterp==0, geometry=3; end
% perform checks on e and w
if nargin>=5
  if any(size(w)==1), w=w(:)'; end
  if size(e,1)~=size(w,2)
    if size(e,1)==size(w,1) && size(w,2)==size(X,1)
    else
      error('w is not compatible in size with e and/or X')
    end
  end
  if any(w<0)
    warnings{end+1}='In g2P: w contains negative elements';
  elseif any(abs(sum(w,2)-1)>size(w,2)*1e-14)
    warnings{end+1}='In g2P: probability weights do not sum to 1';
  end
end
if isnumeric(s)
  if size(s,2)==1
    s={s};
  else
    error('if s is numeric it must be a column vector')
  end
end
% if s is a structure variable several alternative grid geometries are 
% permitted.
if isstruct(s)
  if isfield(s,'type')
    switch s.type
      case 'rectangular'
        geometry=0;
        s=s.grid;
      case 'scatter'
        geometry=1;
      case 'simplex'
        geometry=2;
    end
  else
    geometry=1;  % rectangular
  end
end

if nargin<5 || isempty(e)   % deterministic problem
  P = getPk(feval(g,X),s,geometry,cleanup);
else         % stochastic problem
  try
    if iscell(e) 
      if iscell(w) && length(e)==length(w)
        [e,w]=disccombine(e,w);
      else
        error('if e is a cell array w must be a cell array with the same # of elements')
      end
    end
    nx=size(X,1);
    ne=size(e,1);
    P=getPk(feval(g,repmat(X,ne,1),kron(e,ones(nx,1))),s,geometry,cleanup);
    if any(size(w)==1)
      P = P*kron(w(:),speye(nx));
    else
      w=w';
      P = mxv(P,w(:))*kron(ones(ne,1),speye(nx));
    end
  catch
    if expande
      zz=zeros(nx,1);
    else
      zz=0;
    end
    if size(e,1)==size(w,1) && size(w,2)==size(X,1)
      w=w';
    end
    gval = feval(g,X,e(1+zz,:)); 
    P=getPk(gval,s,geometry,cleanup);
    clear gval;
    ns=size(P,1);
    P = mxv(P,w(:,1));
    K=size(w,2);
    for k=2:K
      gval = feval(g,X,e(k+zz,:));
      Pk=getPk(gval,s,geometry,cleanup);
      clear gval;
      Pk=mxv(Pk,w(:,k));
      P = P + Pk;
      clear Pk;
    end  
  end
end

% eliminate inappropriate values arising from extrapolation
if cleanup==1 && geometry~=1
  try
    P(P<0)=0; 
  catch
    PP=P<0; [ii,jj]=find(PP); 
    for j=1:length(jj), PP(ii(j),jj(j))=0; end
  end
  P=mxv(P,1./sum(P));
else
  if any(any(P<0))
    warnings{end+1}='In g2P: P contains negative elements';
  end
end
if nargout<2
  for i=1:numel(warnings)
    disp(warnings{i})
  end
end


function Pk=getPk(gval,s,geometry,cleanup)
  switch geometry
    case 0
      Pk = rectbas(gval,s,[],cleanup);
    case 1
      Pk = scatterbas(gval,s);
    case 2
      Pk = simplexbas(gval,s.params{:});
    case 3
      Pk = freudenthal(gval,s,cleanup);
  end