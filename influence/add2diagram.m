% add2diagram Adds a variable to an influence diagram
% USAGE
%   D=add2diagram(D,name,type,obs,parents,cpd);
% INPUTS
%   D           : an existing influence diagram structure 
%                  (empty if creating a new diagram)
%   name        : string variable name
%   type        : variable type
%                   's' current state
%                   'a' or 'd' action (decision)
%                   'f' future state
%                   'r' or 'u' reward (utility)
%                   'p' parameter
%                   'c' chance (not one of the above)
%   obs         : 1 if variable is observed [the defualt], 0 if unobserved
%   parents     : cell array of strings listing the variable's parents 
%   cpd         : rv structure (create using rvdef)
%   loc         : 2-vector of numbers on (0,1) specifying location on figure [optional]
%   attachments : 2-column matrix of attachment locations for each parent 
%                   and for own node (empty if variable has no parents)
%                   attachments are in {1,...,8} and specify where on a
%                   node an arc attaches
% OUTPUT
%   D           : updated influence diagram structure

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

function D=add2diagram(D,name,type,obs,parents,cpd,locs,attachments)
if nargin<8
  attachments=[];
if nargin<7
  locs=[];
if nargin<6
  cpd=[];
if nargin<5
  parents={[]};
if nargin<4
  obs=true;
if nargin<3
  type='';
if nargin<2
  error('At least 2 inputs must be passed')
end; 
end; 
end; 
end; 
end; 
end;
end;
if ischar(parents), parents={parents}; end
if isempty(D)  
  n=1;
  D.names={name};
  D.types={type};
  D.obs=obs;
  if isempty(parents)
    D.parents={{}};
  elseif ischar(parents)  % this should never happen because no parents can be defined for the first variable
    D.parents={parents};
  else
    D.parents=parents;
  end
  D.cpds={cpd};
  if nargin>=6
    D.cpdnames={inputname(6)};
  else
    D.cpdnames={[]};
  end
  D.values={};
  D.locs=locs;
  D.attachments=zeros(0,4);
else
  n=length(D.names)+1;
  D.names{n}=name;
  D.types{n}=type;
  D.obs=[D.obs obs];
  if isempty(parents)
    D.parents{n}={};
  elseif ischar(parents)
    D.parents(n)={parents};
  else
    D.parents{n}=parents;
  end
  D.cpds{n}=cpd;
  if nargin>=6
    D.cpdnames{n}=inputname(6);
  else
    D.cpdnames{n}=[];
  end
  D.locs=[D.locs;locs];
  [check,pnum]=ismember(parents,D.names);
  if isempty(attachments) % use default values 5 & 1  
    if pnum>0
      D.attachments=[D.attachments; [pnum' n+zeros(length(pnum),1) 5+zeros(length(pnum),1) ones(length(pnum),1)]];
    end
  else 
    D.attachments=[D.attachments; [pnum' n+zeros(length(pnum),1) attachments]];
  end
end

if ~isempty(cpd)
if isnumeric(cpd)
  if size(cpd,2)==1
    cpd=rvdef('v',[],cpd);
    D.cpds{end}=cpd;
  else
    error('cpd incorrectly defined')
  end
elseif isa(cpd,'function_handle')
  if strcmp(D.types{end},'f')
    ii=find(strcmp(D.names,D.names{end}(1:end-1)));
    if isempty(ii)
      error('values must be specified for future state variables')
    end
    cpd=rvdef('f',cpd,D.values{ii});
  else
    cpd=rvdef('f',cpd);
  end
  D.cpds{end}=cpd;
end
end

if isempty(cpd) 
  D.sizes(n)=0;
  D.values{n}=[];
elseif isstruct(cpd)
  if ~isfield(cpd,'size') || isempty(cpd.size), D.sizes(n)=0; 
  else                                           D.sizes(n)=cpd.size;
  end
  if ~isfield(cpd,'values') || isempty(cpd.values), D.values{n}=[]; 
  else                                              D.sizes(n)=length(cpd.values); D.values{n}=cpd.values;
  end
  if strcmp(cpd.type,'d') && isempty(cpd.parameters) && ~isempty(parents)
    [~,pp]=ismember(parents,D.names);
    D.cpds{end}=getmatchfunc(cpd,D.values(pp));
  end
end

ii=ismember(parents,D.names);
if ~all(ii)
  error('parents must be added to a diagram before children')
end

% perform checks on cpd
if ~isempty(cpd) 
  %rvcheck(cpd,parents)
end
