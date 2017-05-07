% updatediagram Updates (or creates) a diagram from a figure
% USAGE
%   Dout=updatediagram(Din,fig,print,dname);
% INPUTS
%   Din   : an influence diagram structure (if omitted the figure is used)
%   fig   : a figure handle (if omitted the current figure is used)
%   print : 1 to print code containing the location and arc attachment information
%           2 to print code for diagram creation and arc attachment
% OUTPUT
%   Dout  : an influence diagram structure based on information 
%             in Din and fig

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

function DD=updatediagram(D,fig,print,dname,typestring)
if nargin<4
  dname='D';
end
if nargin<2 || isempty(fig)
  fig=gcf;
end
if nargin<1
  D=[];
end
if nargin<3 || isempty(print)
  if isempty(D)
    print=2;
  else
    print=0;
  end
end
figinfo=get(fig,'UserData');
hnode=figinfo.nodes;
if isempty(hnode)
  disp('diagram is empty')
  return
end
d=length(hnode);
% # of variables in the original diagram
if isempty(D)
  k=0;
else
  k=length(D.sizes);
end
g=zeros(1,d);
DD=struct('names',[],'types',[],'obs',[],'parents',[],'cpds',[],...
   'values',[],'sizes',zeros(1,d),'locs',[],'attachments',[]);
DD.values=cell(1,d);
% get node information
for i=1:d
  nodedata=get(hnode(i),'UserData');
  if isempty(D)
    ii=[];
  else
    ii=find(ismember(D.names,get(nodedata.texthandle,'String')));
  end
  if isempty(ii)  % new variable(s) added in figure
    k=k+1;
    ii=k;
    DD.cpds{ii}=[];
    DD.sizes(ii)=0;
  else
    DD.sizes(ii)=D.sizes(ii);
    DD.cpds{ii}=D.cpds{ii};
    DD.cpdnames{ii}=D.cpdnames{ii};
    DD.orders{ii}=D.orders{ii};
  end
  pos=get(nodedata.texthandle,'position');
  DD.locs(ii,:)=pos(1:2);
  DD.names{ii}=get(nodedata.texthandle,'String');
  DD.types{ii}=nodedata.type;
  DD.obs(ii)=nodedata.obs;
  DD.cpdnames{ii}=nodedata.cpd;
  g(i)=ii;
end

% get arc information
for i=1:d
  ii=g(i);
  DD.parents{ii}={};
  nodedata=get(hnode(i),'UserData');
  nn=length(nodedata.inarcs);
  for j=1:nn
    inarcj=nodedata.inarcs(j);
    for k=1:d
      nodedatak=get(hnode(k),'UserData');
      if any(nodedatak.outarcs==inarcj)
        DD.parents{ii}=[DD.parents{ii} DD.names{g(k)}];
        break
      end
    end
  end
end

A=adjacency(DD);
if any(any(tril(A)))
  [order,errflag]=topoorder(A);
  if errflag
    warning('Diagram is not acyclic')
  else
    disp('Variables are not in topological order and will be reordered')
    DD.names=DD.names(order);
    DD.types=DD.types(order);
    DD.obs=DD.obs(order);
    DD.parents=DD.parents(order);
    DD.cpds=DD.cpds(order);
    DD.cpdnames=DD.cpdnames(order);
    DD.values=DD.values(order);
    DD.sizes=DD.sizes(order);
    DD.locs=DD.locs(order,:);
  end
end

% get arc attachment information
ha=findall(fig,'Tag','arc');
na=length(ha);
DD.attachments=zeros(na,4);
for k=1:na
  arcdata=get(ha(k),'UserData');
  nodedata=get(arcdata.tail,'UserData');
  v=get(nodedata.texthandle,'String');
  DD.attachments(k,1)=find(ismember(DD.names,v));
  nodedata=get(arcdata.head,'UserData');
  v=get(nodedata.texthandle,'String');
  DD.attachments(k,2)=find(ismember(DD.names,v));
  DD.attachments(k,3:4)=arcdata.attachments;
end

parents=getparents(DD);
if print
  if print>1
    % print dummy variables for CPDs
    % so code will run even though cpds are not defined
    for i=1:d      
      if ~isempty(DD.cpdnames{i})
        fprintf([DD.cpdnames{i} '=[];\n'])
      end
    end
    k=0;
    fprintf([dname '=[];\n']);
    for i=1:d      
      fprintf([dname '=add2diagram(' dname ',''' DD.names{i} ...
          ''',''' typestring(DD.types{i}) ''','])  
      if DD.obs(i), fprintf('true,{')
      else fprintf('false,{')
      end
      for j=1:length(parents{i})
        if j>1, fprintf(','); end
        fprintf(['''' DD.names{parents{i}(j)} ''''])
      end
      if isempty(DD.cpdnames{i})
        fprintf('},[],')
      else
        fprintf(['},' DD.cpdnames{i} ','])
      end
      fprintf('[%1.4f,%1.4f]',DD.locs(i,1),DD.locs(i,2));
      
      if isempty(DD.attachments)
        fprintf(');\n')
      else
        fprintf(',[')
        ii=find(DD.attachments(:,2)==i);
        for j=1:length(parents{i})
          k=find(DD.attachments(ii,1)==parents{i}(j)); 
          if j==1
            fprintf('%1i %1i',  DD.attachments(ii(k),3), DD.attachments(ii(k),4))
          else
            fprintf(';%1i %1i', DD.attachments(ii(k),3), DD.attachments(ii(k),4))
          end
        end
        fprintf(']);\n')
      end
    end
  end
  
  % location and attachment info now in node line
  if 0
    disp([dname '.locs=[ ...'])
    for j=1:2
      for i=1:size(DD.locs,1)
        if i>1, fprintf(' '); end
        fprintf('%5.3f',DD.locs(i,j));
      end
      if j<2, fprintf(';\n')
      else    fprintf(']'';\n')
      end
    end
    if ~isempty(DD.attachments)
      disp([dname '.attachments=[ ...'])
      for j=1:4
        for i=1:size(DD.attachments,1)
          if i>1, fprintf(' '); end
          fprintf('%2i',DD.attachments(i,j));
        end
        if j<4, fprintf(';\n')
        else    fprintf(']'';\n')
        end
      end
    end
  end
  
  displayfigoptions

end

 function displayfigoptions
%   name            : name appearing at the top of the figure 
%   figpos          : position of figure
%   fontsize        : size of label font [0.035] in normalized units
%   fontname        : a valid font name ['Rockwell'] for all text
%                       type 'listfonts' for available fonts
%   backgroundcolor : background color of plot [light gray - [.85 .85 .85]]
%   nodecolor       : default node color [ light blue - [0 .95 .9]]

figinfo=get(gcf,'UserData');
fprintf('foptions=struct(...\n')
fprintf('''name'',''')
fprintf('%1c',figinfo.name);
fprintf(''',...\n')
fprintf('''figpos'',[%5.3f %5.3f %5.3f %5.3f],...\n',get(gcf,'position'));
fprintf('''fontsize'',%1.3f,...\n',figinfo.fontsize);
fprintf('''fontname'',''')
fprintf('%1c',figinfo.fontname);
fprintf(''',...\n')
fprintf('''backgroundcolor'',[%5.3f %5.3f %5.3f],...\n',figinfo.backgroundcolor);
fprintf('''nodecolor'',[%5.3f %5.3f %5.3f]);\n',figinfo.nodecolor(1,:));

return

fprintf('foptions=struct(...\n') %#ok<UNRCH>
if length(figinfo.name)>1
  fprintf('''name'',''')
  fprintf('%1c',figinfo.name);
  fprintf(''',...\n')
end

fprintf('''figpos'',[%5.3f %5.3f %5.3f %5.3f],...\n',get(gcf,'position'));

if figinfo.fontsize ~= 0.035
  fprintf('''fontsize'',%1.3f,...\n',figinfo.fontsize);
end
if ~strcmp(figinfo.fontname,'Rockwell')
  fprintf('''fontname'',''')
  fprintf('%1c',figinfo.fontname);
  fprintf(''',...\n')
end
if any(figinfo.backgroundcolor~=[.85 .85 .85])
  fprintf('''backgroundcolor'',[%5.3f %5.3f %5.3f],...\n',figinfo.backgroundcolor);
end
if any(figinfo.nodecolor~=[0 .95 .9])
  fprintf('''nodecolor'',[%5.3f %5.3f %5.3f],...\n',figinfo.nodecolor);
end
fprintf(');\n')

