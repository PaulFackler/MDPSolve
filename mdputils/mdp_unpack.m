% mdp_unpack Gets information from an MDPSOLVE model variable, converts to a standard form and performs checks

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

function [errors,warnings,R,P,delta,ns,nx,Ix,Iexpand,colstoch,EV,T,nstage,nrep,Xindexed,expandP] ...
             = mdp_unpack(model,debug)        
  errors={};
  warnings={};
  if ~isstruct(model)
    if debug, error(' '); end
    errors{end+1}={65};
    return
  end
  
  celltype=true; numtype=false; % used when calling getmodelfield
  % checks if a vector contains all positive integers; note: empty and inf return true
  positiveinteger = @(x) isnumeric(x) && all((x(:)>=1) & (x(:)==floor(x(:))));
  % check if x is a non-empty vector or matrix of numbers; can't be logical array
  ismatrix = @(x) ~isempty(x) && isnumeric(x) && ndims(x)<=2;
  
  %%%% horizon or T field
  if isfield(model,'horizon')
    T = model.horizon;
    if isfield(model,'T')
      if debug, error(' '); end
      errors{end+1}=71;
    end
  elseif isfield(model,'T')
    T = model.T;
  else
    T=inf;
  end
  if ~positiveinteger(T) || numel(T)~=1
    if debug, error(' '); end
    errors{end+1}=72;
  end
 
  % nstage field
  if isfield(model,'nstage')
    nstage=model.nstage;
    if ~positiveinteger(nstage) || numel(nstage)~=1
      if debug, error(' '); end
      errors{end+1}=81;
      nstage=[];
    end
  else
    nstage=[];
  end
  
   %%%% reward or R field
  [R,nstage]=getmodelfield(model,'R',nstage,celltype);
  if isempty(R)
    [R,nstage]=getmodelfield(model,'reward',nstage,celltype);
  else
    if isfield(model,'reward')
      if debug, error(' '); end
      errors{end+1}=73;
    end
  end
  for i=1:numel(R)
    if ~isa(R{i},'function_handle') && ~ismatrix(R{i})
      if debug, error(' '); end
      errors{end+1}=74; %#ok<AGROW>
    end
  end
  
  %%%%% discount rate
  if ~isfield(model','d') && ~isfield(model','discount')
      if debug, error(' '); end
      errors{end+1}=8;
  end
  [delta,nstage]=getmodelfield(model,'discount',nstage,celltype);
  if isempty(delta)
    [delta,nstage]=getmodelfield(model,'d',nstage,celltype);
  else
    if isfield(model,'d')
      if debug, error(' '); end
      errors{end+1}=75;
    end
  end
  if ~isempty(delta)
    for i=1:numel(delta)
      if ~isa(delta{i},'function_handle') && ~ismatrix(delta{i})
        if debug, error(' '); end
        errors{end+1}=76;  %#ok<AGROW>
      elseif all(size(delta{i})>1)
        if debug, error(' '); end
        errors{end+1}=76;  %#ok<AGROW>
      end
    end
  end
  
  %%%% transition probability matrices
  [P,nstage]=getmodelfield(model,'transprob',nstage,celltype);
  if isempty(P)
    [P,nstage]=getmodelfield(model,'P',nstage,celltype);
  else
    if isfield(model,'P')
      if debug, error(' '); end
      errors{end+1}=77;
    end
  end
  
  for i=1:numel(P)
    if ~isa(P{i},'function_handle') && ~ismatrix(P{i})
      if debug, error(' '); end
      errors{end+1}=78; %#ok<AGROW>
    elseif ismatrix(P{i}) 
      if any(size(P{i})==1)
        if debug, error(' '); end
        errors{end+1}=78; %#ok<AGROW>
      elseif any(any(isnan(P{i})))
        if debug, error(' '); end
        errors{end+1}=26; %#ok<AGROW>
      end
    end
  end
  
  % determine number of state values (ns)
  [ns,nstage]       = getmodelfield(model,'ns',nstage,numtype);
    if ~isempty(ns) && any(ns<1 | ns~=floor(ns))
      if debug, error(' '); end
      errors{end+1}=79;
    end
  [nx,nstage]       = getmodelfield(model,'nx',nstage,numtype);
    if ~isempty(nx) && any(nx<1 | nx~=floor(nx))
      if debug, error(' '); end
      errors{end+1}=79;
    end
  [X,nstage]        = getmodelfield(model,'X',nstage,celltype);
  [svars,nstage]    = getmodelfield(model,'svars',nstage,celltype);
  [Ix,nstage]       = getmodelfield(model,'Ix',nstage,celltype);
  [nrep,nstage]     = getmodelfield(model,'nrep',nstage,numtype);
    if isempty(nrep), nrep=1; end
    if ~positiveinteger(nrep)
      if debug, error(' '); end
      errors{end+1}=80;
    end
  % P is a function to compute E[V(S+)|X]
  [EV,nstage]       = getmodelfield(model,'EV',nstage,numtype);
    if isempty(EV)
      EV=false;
    else
      EV=logical(EV);
    end
  % specifies if probability matrices have ns rows or ns columns
  [colstoch,nstage] = getmodelfield(model,'colstoch',nstage,numtype);
  [Iexpand,nstage]  = getmodelfield(model,'Iexpand',nstage,celltype);
      
  % after going through all the valid fields we should know nstage
  % if it is not yet determined it must be 1
  if isempty(nstage)
    nstage=1; 
  elseif ~positiveinteger(nstage) % can't be true can it?
    if debug, error(' '); end
    errors{end+1}=81;
  end
 
  
  % This is no longer needed ???
  % check all of the sizes
  err=checksize(P,nstage);
  if err, if debug, error(' '); end; errors{end+1}={82,'P'}; end
  err=checksize(R,nstage);
  if err, if debug, error(' '); end; errors{end+1}={82,'R'}; end
  err=checksize(delta,nstage);
  if err, if debug, error(' '); end; errors{end+1}={82,'d'}; end
  err=checksize(ns,nstage);
  if err, if debug, error(' '); end; errors{end+1}={82,'ns'}; end
  err=checksize(Ix,nstage);
  if err, if debug, error(' '); end; errors{end+1}={82,'Ix'}; end
  err=checksize(Iexpand,nstage);
  if err, if debug, error(' '); end; errors{end+1}={82,'Iexpand'}; end
  err=checksize(nrep,nstage);
  if err, if debug, error(' '); end; errors{end+1}={82,'nrep'}; end
 
  % determine nx from Ix, Iexpand or R
  if ~isempty(nx), nnx=nx;
  else             nnx=[];
  end
    nn=max([numel(Ix),numel(Iexpand),numel(R)]);
    % if Ix, Iexpand and R are singletons
    if nn==1
      nx=0;
      if ~isempty(Ix)
        nx=numel(Ix{1});
      end
      if nx==0 && ~isempty(Iexpand)
        nx=numel(Iexpand{1});
      end
      if nx==0
        Ri=R{1};
        if ~isnumeric(Ri), Ri=Ri(); end  % Ri is a function handle
        nx=numel(Ri);
        clear Ri
      end
    % otherwise at least one is stage specific
    else
      nx=zeros(1,nn);
      for i=nn:-1:1
        if ~isempty(Ix)
          Ixi=Ix{min(i,numel(Ix))};
          if ~isempty(Ixi)
            nx(i)=numel(Ixi);
          end
          clear Ixi
        end
        if nx(i)==0 && ~isempty(Iexpand)
          Iexpandi=Iexpand{min(i,numel(Iexpand))};
          if ~isempty(Iexpandi)
            nx(i)=numel(Iexpandi);
          end
          clear Iexpandi
        end
        if nx(i)==0
          Ri=R{min(i,numel(R))};
          if ~isnumeric(Ri), Ri=Ri(); end  % Ri is a function handle
          nx(i)=numel(Ri);
          clear Ri
        end
      end
    end
  if ~isempty(nnx) && any(nx~=nnx) 
    if debug, error(' '); end; 
    errors{end+1}={25};
    clear nxx
  end
  if ~positiveinteger(nx)
    if debug, error(' '); end; 
    errors{end+1}=62;
  end
  if numel(nx)==1, nx=nx+zeros(1,nstage); end
  
    
  % determine ns from Ix, R or X and svars
  if ~isempty(ns), nns=ns;
  else             nns=[];
  end
    nn=max([numel(Ix),numel(X),numel(svars),numel(R)]);
    % if Ix, R, X and svars are singletons
    if nn==1
      ns=0;
      Xindexed=true;
      if ~isempty(Ix)
        ns=max(Ix{1});
      end
      if ns==0
        Ri=R{1};
        if ~isnumeric(Ri), Ri=Ri(); end    % Ri is a function handle
        if size(Ri,1)<nx(1)
          ns=size(Ri,1); 
          Xindexed=false;
        end
        clear Ri
      end
      if ns==0 
        if ~isempty(X) && ~isempty(svars)
          Ix{1}=getI(X{1},svars{1});
          ns=max(Ix{1});
        else
          if debug, error(' '); end; 
          errors{end+1}={66};
        end
      end
    % otherwise at least one is stage specific
    else
      ns=zeros(1,nn);
      Xindexed=repmat(true,1,nn);
      for i=nn:-1:1
        if ~isempty(Ix)
          Ixi=Ix{min(i,numel(Ix))};
          if ~isempty(Ixi)
            ns(i)=max(Ixi);
            clear Ixi;
          end
        end
        if ns(i)==0 
          Ri=R{min(i,numel(R))};
          if ~isnumeric(Ri), Ri=Ri(); end   % Ri is a function handle
          if size(Ri,1)<nx(1)
            ns(i)=size(Ri,1); 
            Xindexed(i)=false;
          end
          clear Ri
        end
        if ns(i)==0 
          if ~isempty(X) && ~isempty(svars) ...
          && ~isempty(X{min(i,numel(X))}) && ~isempty(svars{min(i,numel(svars))})
            Ix{i}=getI(X{min(i,numel(X))},svars{min(i,numel(svars))});
            ns(i)=max(Ix{i});
          else
            if debug, error(' '); end; 
            errors{end+1}={67,i}; %#ok<AGROW>
          end
        end
      end
    end
    
  if ~isempty(nns) && any(ns~=nns)  
    if debug, error(' '); end; 
    errors{end+1}={25};
    clear nns
  end
  if ~positiveinteger(ns)
    if debug, error(' '); end; 
     errors{end+1}=63;
  end
  if numel(ns)==1, ns=ns+zeros(1,nstage); end
 
  % perform size checks on Iexpand
  if isempty(Iexpand)
    expandP=false;
  else 
    expandP=repmat(false,1,numel(Iexpand));
    for i=numel(Iexpand):-1:1
      if ~isempty(Iexpand{i})
        if numel(Iexpand{i})~=nx(i)
          if debug, error(' '); end; 
          errors{end+1}={91,i,nx(i)}; %#ok<AGROW>
        end
        expandP(i)=true;
      end
    end
  end
  if numel(expandP)==1; 
    expandP=repmat(expandP,1,nstage);
  end
  
  if numel(Xindexed)==1; Xindexed=repmat(Xindexed,1,nstage); end
    
  % convert to numerical empty
  if isempty(Ix),      Ix=[];      end
  if isempty(Iexpand), Iexpand=[]; end
  
  % check if P contains function handles when EV is true
  for i=max(numel(P),numel(EV)):-1:1
    if EV(min(i,numel(EV))) && ~isa(P{min(i,numel(P))},'function_handle')
      if debug, error(' '); end; 
      errors{end+1}={92,i}; %#ok<AGROW>
    end
  end
  if numel(EV)==1
    EV=repmat(EV,1,nstage);
  end
  
  % determine colstoch 
   if isempty(colstoch)
    if numel(P)==1
      if any(ns~=ns(1)) || any(nx~=nx(1))
        if debug, error(' '); end; 
        errors{end+1}={};
      else
        if EV(1)
          colstoch=true;  % arbitrary as it is not used
        else
          if isnumeric(P{1}), Pi=P{1};
          else                Pi=P{1}();
          end
          if any(any(isnan(Pi)))
            if debug, error(' '); end; 
            errors{end+1}={26};
            code=-2;
          elseif expandP(1)
            code=checkP(Pi,ns(1),max(Iexpand{1}(:)));
          else
            code=checkP(Pi,ns(1),nx(1));
          end
          clear Pi
          switch code
            case -2, if debug, error(' '); end; errors{end+1}={21};
            case -1, if debug, error(' '); end; errors{end+1}={23};
            case  0, colstoch=false;
            case  1, colstoch=true;
          end
        end
      end
      colstoch=repmat(colstoch,1,nstage);
    else
      for i=numel(P):-1:1
        if EV(i)
          colstoch(i)=true;  % arbitrary as it is not used
        else
          if isnumeric(P{i}), Pi=P{i};
          else                Pi=P{i}();
          end
          if any(any(isnan(Pi)))
            if debug, error(' '); end; 
            errors{end+1}={26};
          elseif expandP(i)
            if i>1
              code=checkP(Pi,ns(i-1),max(Iexpand{min(i,numel(Iexpand))}(:)));
            else
              code=checkP(Pi,ns(end),max(Iexpand{min(i,numel(Iexpand))}(:)));
            end
          else
            if i>1
              code=checkP(Pi,ns(i-1),nx(i));
            else
              code=checkP(Pi,ns(end),nx(i));
            end
          end
          clear Pi
          switch code
            case -2, if debug, error(' '); end; errors{end+1}={22,i}; %#ok<AGROW>
            case -1, if debug, error(' '); end; errors{end+1}={24,i}; %#ok<AGROW>
            case  0, colstoch(i)=false;
            case  1, colstoch(i)=true;
          end
        end
      end
    end    
  else
    if numel(colstoch)==1, colstoch=repmat(colstoch,1,nstage); end
  end
  
  % check that delta is scalar or vector of size ns
  for i=1:numel(delta)
    di=delta{i};
    if i<nstage, nsi=ns(i+1); else nsi=ns(1); end
    if ~isnumeric(di), di=di(); end
    if numel(di)>1 && numel(di)~=nsi
      if debug, error(' '); end; 
      errors{1,end+1}={7,i}; %#ok<AGROW>
    end
  end
  if numel(delta)==1 && numel(di)>1 && any(ns~=ns(1))
    if debug, error(' '); end; 
    errors{1,end+1}={7,1};
  end
  clear di
  
  
  % gets a field from model and determines if it either specifies
  % the number of stages or is consistent with previous specification
  % of nstage. vartype=true for cell array outputs and false for numerical array
  % outputs
  function [outvar,nstage]=getmodelfield(model,fname,nstage,vartype)
  if isfield(model,fname)
    outvar=model.(fname);
    if iscell(outvar) && numel(outvar)>1
      if isempty(nstage)
        nstage=numel(outvar);
      else
        if nstage~=numel(outvar)
          if debug, error(' '); end; 
          errors{end+1}={93,fname,nstage};
          if vartype, outvar={[]}; else outvar=[]; end
          return
        end
      end
    end
  else
    if vartype, outvar={}; else outvar=[]; end
    return
  end
  if ~iscell(outvar) && vartype==1
    outvar={outvar};
  end
  % convert cell array to vector
  if vartype==0
    if iscell(outvar)
      temp=outvar;
      outvar=zeros(1,numel(temp));
      for ii=numel(temp):-1:1
        if numel(temp{i})~=1
          outvar=[];
          return
        else
          outvar(ii)=temp{ii};
        end
      end
    end
  end
  end
end
  
% checks to see if a variable is either scalar or has nstage elements
function err=checksize(x,nstage)
  N=numel(x);
  if N>1 && N~=nstage
    err=true;
  else
    err=false;
  end
end
      
  