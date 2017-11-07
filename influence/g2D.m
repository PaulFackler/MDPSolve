% g2D Converts a model with transition functions into a diagrammatic model
% USAGE
%   D=g2D(g,slist,x,e,xelist,names,options);
% INPUTS  
%   g       : ds+1-element cell array of function handles for state transition
%               and utility functions
%   slist   : ds-element vectors of state variable indices in x 
%   x       : dx-element cell array of state/action value vectors
%   e       : de-element cell array of rv structures [optional]
%   xelist  : ds+1-element cell array of lists of variable indices 
%               for state transitios and utility functions
%               element i should contain a list specifying which variables
%               in x and e should be passed to g{i}. Specify x variables
%               with positive integers and e variables with negative integers.
%               Variables must be listed in the order defined by g{i}.
%   names   : dx+de cell array with variables names of the x variables first, 
%               then the e variables
%   options : a structure variable to control the procedure
%                'print' =1 to print out the diagram code to cut & paste 
% OUTPUTS
%   D       : a stucture variable specifying a diagrammatic model
%
% The functions in g should accept a k-row column vectors of states, action and 
%   noise variables and return a k-row column vector of next period state values


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

function D=g2D(g,slist,x,e,xelist,names,options)
print=0;
order=[];
if nargin<7, options=[]; end
if isnumeric(options)
  print=options;
else
  getopts(options, ...
         'print',       0, ...
         'order',       []); 
end
g=inputcheck(g,'g','function_handle');
x=inputcheck(x,'x','numeric');
e=inputcheck(e,'e','struct');
xelist=inputcheck(xelist,'xlist','numeric');
s=x(slist);

ds=length(s);
dx=length(x);
de=length(e);

avars=~ismember(1:dx,slist);
svars=find(~avars);
avars=find(avars);

xlist=cell(1,ds);
elist=cell(1,ds);
for i=1:ds
  if length(xelist{i})~=nargin(g{i})
    error(['xelist has a different number of inputs than g for variable ' num2str(i)])
  end
  elist{i}=-xelist{i}(xelist{i}<0);
  if any(elist{i}>de)
    error(['xelist{' num2str(i) '} contains a value greater than the # of noise variables'])
  end
  xlist{i}=xelist{i}(xelist{i}>0);
  if any(xlist{i}>dx)
    error(['xelist{' num2str(i) '} contains a value greater than the # of X variables'])
  end
end
if isempty(names)
  for i=1:ds,    names{svars(i)}=['S' num2str(i)];  end
  for i=1:dx-ds, names{avars(i)}=['A' num2str(i)];  end
  for i=1:de,    names{dx+i}    =['e' num2str(i)];   end
end
cpdnames=cell(1,dx+de);
for i=1:ds,    cpdnames{i}        =['x{' num2str(svars(i)) '}'];    end
for i=1:dx-ds, cpdnames{ds+i}     =['x{' num2str(avars(i)) '}']; end
for i=1:de,    cpdnames{dx+i}     =['e{' num2str(i) '}'];    end
for i=1:ds,    cpdnames{dx+de+i}  =['g{' num2str(i) '}'];    end
D=[];
for i=1:ds
  D=add2diagram(D,names{i},'s',1,{},s{i});
end
for i=1:dx-ds
  D=add2diagram(D,names{ds+i},'a',1,{},x{ds+i});
end
for i=1:de
  D=add2diagram(D,names{dx+i},'c',1,{},e{i});
end
for i=1:ds
  di=length(xelist{i});
  parentsi=cell(1,di);
  for j=1:di
    if xelist{i}(j)>0, parentsi{j}=names{xelist{i}(j)};
    else               parentsi{j}=names{dx-xelist{i}(j)};
    end
  end
  D=add2diagram(D,[names{i} '+'],'f',1,parentsi,g{i});
end
if length(g)>ds
  i=ds+1;
  di=length(xelist{i});
  parentsi=cell(1,di);
  for j=1:di
    if xelist{i}(j)>0, parentsi{j}=names{xelist{i}(j)};
    else               parentsi{j}=names{dx-xelist{i}(j)};
    end
  end
  D=add2diagram(D,['U'],'u',1,parentsi,g{i});
  cpdnames{dx+de+ds+1}   =['g{' num2str(i) '}'];
end
D.cpdnames=cpdnames;
if print
  printD(D,1);
end

function x=inputcheck(x,xname,class)
if ~iscell(x)
  if isa(x,'class')
    x={x};
  else
    error([xname ' must be a single element or a cell array of elements of class ' class])
  end
end