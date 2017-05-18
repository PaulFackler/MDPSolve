% drawdiagram Creates a plot of a diagram and allows interactive editing
% USAGE
%   drawdiagram(D,options)
% INPUTS
%   D       : an influence diagram structure ([] for a blank diagram)
%   options : a structure variable to control figure properties
%
% Actions are displayed as rectangles
% The utility variable is displayed as a diamond
% All other variables are displayed as ellipses
% Variables defined by a deterministic function have thick boundaries
%
% The figure can be manipulated interactively.
%   Right click to select or deselect a node or arc
%   Right click a second node to add an arc from the selected node 
%   Delete key to remove a selected node or arc
%   Left click and hold to move a node or arc
%   Arcs can only be moved to one of the 8 attachment points on a node; 
%     a moved arc will snap to the nearest attachment point
%   Left click on a blank area to add a node there; a dialog box will open
%   F1 Help screen
%   F2 change node color
%   F3 change background color
%   F4 change font name
%   F5 toggle mean and tool bar
%   F6 output graph information to screen
%   F7 toggle confirmation boxes
%   Arrows, PgUp, PgDn, Home and End change figure size and placement
%
% Options: 
%   name            : name appearing at the top of the figure 
%   fontsize        : size of label font [0.035] in normalized units
%   fontname        : a valid font name ['Rockwell'] for all text
%                       type 'listfonts' for available fonts
%   addnumbers      : add numbers to labels [false]
%   backgroundcolor : background color of plot [light gary - [.85 .85 .85]]
%   nodecolor       : default node color [ light blue - [0 .95 .9]]
%   obsfactor       : factor on [0,1] to reduce saturation for 
%                       non-observed variables [0.15]
%   nodeedgewidth   : width of the line defining the nodes [1]
%   arcflag         : 0 arcs are connected at closest points
%                     1 arcs go from rightmost point to leftmost point 
%                     2 arcs go from bottom point to top point [0]
%   arcwidth        : width of the arcs [2]
%   archeadsize     : size of arc arrow heads [1]
%   confirmationboxes : if true user must confirm certain choices [false]
%                          use F7 to toggle this
%   undirected        : true for undirected graphs (superceeds dag) [false]
%   dag               : true for DAG (Directed Acyclic Graph) [true]
%   generic           : true for non-MDPSolve graph [false]
%   nodefunc          : handle to a node edit box  '
%   shapes            : k-vector on {1,...,7} mapping types to shapes 
%   typestring        : k-vector of letter designations for types
%   types             : q-vector (q<=k) of integers 
  
%
% All units are normalized to the figure window
% All colors are RGB triplets (see MATLAB color documentation)
% Each node has 8 possible attachment points for arcs. 
% When arcflag=1 attachments are made between the nearest points between
%   the parent and child nodes. When moving an arc the attachment is made
%   to the nearest attachment point on the appropriate node.
%
% By default this procedure turns off Matlab's Menu Bar and Tool Bar features
% To restore these use
%   set(fig,'MenuBar','figure')
% or
%   set(fig,'ToolBar','figure')
% where fig is the figure number
%
% Use updatediagram (or F6) to extract the information from the figure

% Technical notes:
% Information about general options is stored in fig.UserData
% Each variable (node) has a structure associated with it 
%   that is stored in handle.UserData and contains fields:
%     texthandle : handle to the text object
%     type       : variable type (positive integer)
%                     default is s=1, a,d=2, f=3, u,r=4, c=5, p=6
%     obs        : true/false indicator if variable is observed
%     cpd        : string that defines the Conditional Probability Distribution
%     inarcs     : vector of arc handles entering node
%     outarcs    : vector of arc handles exiting node
% 
% Each arc has a structure associated with it 
%   that is stored in handle.UserData and contains fields:
%     head : handle of the node associated with the arc head (child)
%     tail : handle of the node associated with the arc tail (parent)
%     attachments : 2-vector with the attachment points for the 
%                     tail and head nodes
%
% In addition to fields for the general options fig.UserData contains 
%   the following fields:
%     nodes : a list of all node handles
%     arcs  : a list of all arc handles

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

function drawdiagram(D,options)
cf=gcf;
set(cf,'Visible','off')
clf('reset');
set(cf,'units','normalized')

figinfo.name='';
figinfo.figpos=get(gcf,'position');
figinfo.fontsize = 0.035;
figinfo.fontname = 'Rockwell';
figinfo.addnumbers=false;
figinfo.backgroundcolor=[.85 .85 .85];
figinfo.obsfactor = 0.15;
figinfo.nodeedgewidth=1;
figinfo.arcflag=0;
figinfo.arcwidth=2;
figinfo.archeadsize=0.75;
figinfo.boxsize=1;
figinfo.confirmationboxes=false;
figinfo.undirected=false;
figinfo.dag=true;
figinfo.generic = false;
figinfo.nodeedit = @getnewvariable;
figinfo.typestring='sadfurcp';    % valid letters to designate types
figinfo.types=[1 2 2 3 4 4 5 6];  % mapping from letter types to type #
figinfo.shapes=[2 1 2 4 5 3];     % mapping from type # to shape #
if nargin>=2 && ~isempty(options)
  if isfield(options,'figpos'),            figinfo.figpos=options.figpos;                       end
  if isfield(options,'name'),              figinfo.name=options.name;                           end
  if isfield(options,'fontsize'),          figinfo.fontsize=options.fontsize;                   end
  if isfield(options,'fontname'),          figinfo.fontname=options.fontname;                   end
  if isfield(options,'addnumbers'),        figinfo.addnumbers=options.addnumbers;               end
  if isfield(options,'backgroundcolor'),   figinfo.backgroundcolor=options.backgroundcolor;     end
  if isfield(options,'nodecolor'),         figinfo.nodecolor=options.nodecolor;                 end
  if isfield(options,'statecolor'),        figinfo.statecolor=options.statecolor;               end
  if isfield(options,'actioncolor'),       figinfo.actioncolor=options.actioncolor;             end
  if isfield(options,'utilitycolor'),      figinfo.utilitycolor=options.utilitycolor;           end
  if isfield(options,'chancecolor'),       figinfo.chancecolor=options.chancecolor;             end
  if isfield(options,'paramcolor'),        figinfo.paramcolor=options.paramcolor;              end
  if isfield(options,'obsfactor'),         figinfo.obsfactor=options.obsfactor;                end
  if isfield(options,'nodeedgewidth'),     figinfo.nodeedgewidth=options.nodeedgewidth;         end
  if isfield(options,'arcflag'),           figinfo.arcflag=options.arcflag;                     end
  if isfield(options,'arcwidth'),          figinfo.arcwidth=options.arcwidth;                   end
  if isfield(options,'archeadsize'),       figinfo.archeadsize=options.archeadsize;             end
  if isfield(options,'confirmationboxes'), figinfo.confirmationboxes=options.confirmationboxes; end
  if isfield(options,'undirected'),        figinfo.undirected=options.undirected;               end
  if isfield(options,'dag'),               figinfo.dag=options.dag;                             end
  if isfield(options,'generic'),           figinfo.generic=options.generic;                     end
  if isfield(options,'nodeedit'),          figinfo.nodeedit=options.nodeedit;                   end
  if isfield(options,'shapes'),            figinfo.shapes=options.shapes;                       end
  if isfield(options,'typestring'),        figinfo.typestring=options.typestring;               end
  if isfield(options,'types'),             figinfo.types=options.types;                         end
  if isfield(options,'modecolor'),         figinfo.nodecolor=options.nodecolor;                 end
  %if isfield(options,''),  figinfo.=options.;   end
end
if figinfo.undirected, figinfo.dag=false; figinfo.archeadsize=0; end
figinfo.addnode=@addnode;
figinfo.unselect=@unselect;

if ~isfield(figinfo,'nodecolor')
  nodecolor   = repmat([0 .95 .9],length(figinfo.shapes),1);
  if isfield(figinfo,'statecolor'),   
    nodecolor(1,:)=figinfo.statecolor;
    nodecolor(3,:)=figinfo.statecolor;
  end
  if isfield(figinfo,'chancecolor'),   
    nodecolor(4,:)=figinfo.chancecolor;
  end
  if isfield(figinfo,'paramcolor'),   
    nodecolor(6,:)=figinfo.paramcolor;
  end
  if isfield(figinfo,'actioncolor'),   
    nodecolor(2,:)=figinfo.actioncolor;
  end
  if isfield(figinfo,'utilitycolor'),   
    nodecolor(4,:)=figinfo.utilitycolor;
  end
  figinfo.nodecolor=nodecolor;
end

set(gcf,'color',figinfo.backgroundcolor,...
    'Position',figinfo.figpos,...
    'KeyPressFcn',@keypressfcn,'ButtonDownFcn',@mouseclick,...
    'ResizeFcn',@redraw,...
    'Name',figinfo.name,'NumberTitle','on','MenuBar','none','ToolBar','none',...
     'PaperPositionMode','auto',...
    'Interruptible','off');

axes('Units','normalized','position',[0 0 1 1],...
     'xtick',[],'ytick',[],...
     'color',figinfo.backgroundcolor,...
     'xcolor',figinfo.backgroundcolor,'ycolor',figinfo.backgroundcolor,...
     'FontUnit','normalized','FontSize',figinfo.fontsize,...
     'xlim',[0 1],'ylim',[0 1],...
     'ButtonDownFcn',@mouseclick,...
     'Visible','off');
set(cf,'UserData',figinfo);
drawnow

if nargin==0 || isempty(D)
  figinfo.nodes=[];
  figinfo.arcs=[];
else
  % set default locations
  if ~isfield(D,'locs') || isempty(D.locs)
    d=length(D.names);
    locs=zeros(d,2);
    nrows=ceil(sqrt(d));
    locs(:,1)=ceil((1:d)/nrows)';
    locs(:,2)=(nrows+1)-((1:d)'-(locs(:,1)-1)*nrows);
    D.locs=locs/(nrows+1);
    clear d locs nrows
    D.attachments=[];
  end
  d=length(D.names);
  figinfo.nodes=zeros(d,1);
  % draw the nodes
  for i=1:d
    type=D.types{i};
    if ischar(type)
      type=figinfo.types(figinfo.typestring==type);
    end
    if isfield(D,'cpdnames')
      figinfo.nodes(i)=addnode(D.locs(i,:),D.names{i},type,D.obs(i),D.cpdnames{i});
    elseif isfield(D,'obs')
      figinfo.nodes(i)=addnode(D.locs(i,:),D.names{i},type,D.obs(i),[]);
    else
      figinfo.nodes(i)=addnode(D.locs(i,:),D.names{i},type,true,[]);
    end
  end
  set(cf,'UserData',figinfo);

  % draw the arcs
  parents=getparents(D);
  na=length([parents{:}]);  % # of arcs
  figinfo.arcs=zeros(na,1);  % not actually used
  % check if attachment information is correct
  attachok=1;
  if isfield(D,'attachments') && ~isempty(D.attachments)
    for k=1:na
        attachments=D.attachments(k,:);
        if ~any(parents{attachments(2)}==attachments(1))
          attachok=0; break; 
        end
    end
  end
  
  if isfield(D,'attachments') && ~isempty(D.attachments) && attachok
    for k=1:na
        attachments=D.attachments(k,:);
        hTail=figinfo.nodes(attachments(1));
        hHead=figinfo.nodes(attachments(2));
        figinfo.arcs(k)=addarc(hTail,hHead,attachments(3:4));
    end
  else
    k=0;
    for i=1:d
      if ~isempty(parents{i}) && D.types{i}~='n' 
        [inset,locs]=ismember(parents{i},1:d);
        if ~all(inset)
          error(['parent list for ' D.names{i} ' contains names not in list'])
        end
        hHead=figinfo.nodes(i);
        for j=1:length(locs)
          k=k+1;
          hTail=figinfo.nodes(locs(j));
          figinfo.arcs(k)=addarc(hTail,hHead);
        end
      end
    end
  end
end
set(cf,'UserData',figinfo);

drawnow
set(cf,'Visible','on');

  
function h2=addnode(loc,name,type,obs,cpd,h2)
  cf=gcf;
  figinfo=get(cf,'UserData');
  nodecolor=figinfo.nodecolor;
  if size(nodecolor,1)>1, nodecolor=nodecolor(type,:); end
  if ~obs, 
    nodecolor=1-(1-nodecolor)*figinfo.obsfactor; 
  end
  
  % define text object
  %loc(1)=ceil(loc(1)*512)/512; % helps prevent leftward drift due to rounding
  h1 = text(loc(1),loc(2),name);
  set(h1,'BackgroundColor',nodecolor,...
      'HorizontalAlignment','center','VerticalAlignment','middle',...
      'FontUnits','normalized','FontSize',figinfo.fontsize,...
      'FontName',figinfo.fontname,...
      'LineWidth',figinfo.nodeedgewidth,...
      'Units','normalized','margin',1);
  %set(h1,'Units','pixels'); set(h1,'Units','data');
  ex=get(h1,'extent');
  if length(name)<3, ex(1)=ex(1)-ex(3)/2; ex(3)=ex(3)*2; end
  ex=ex+[-0.15*ex(3) -0.15*ex(4) 0.3*ex(3) 0.3*ex(4)];
  set(h1,'ButtonDownFcn',@mouseclick);
  [x,y]=getshape(ex,figinfo.shapes(type));
  if nargin<6 || isempty(h2)
    h2=patch(x,y,nodecolor);
    nodeinfo.inarcs=[];
    nodeinfo.outarcs=[];
    set(h2,'ButtonDownFcn',@mouseclick,'tag','node');
  else
    nodeinfo=get(h2,'UserData');
    delete(nodeinfo.texthandle);
    set(h2,'XData',x,'Ydata',y,'FaceColor',nodecolor);
  end
  %if ~isempty(cpd) && cpd==2, set(h2,'LineWidth',figinfo.nodeedgewidth*2); end
  nodeinfo.texthandle=h1;
  nodeinfo.type=type;
  nodeinfo.obs=obs;
  nodeinfo.cpd=cpd;
  set(h2,'UserData',nodeinfo)
  set(h1,'Userdata',h2,'Tag','nodetext')  % store patch handle with text object
  uistack(h1,'bottom')
  uistack(h2,'bottom')

function   [x,y]=getshape(ex,shape)
switch shape
case 1  % rectangle
x=[ex(1)         ex(1) ex(1)+ex(3)/2  ex(1)+ex(3) ex(1)+ex(3)   ex(1)+ex(3)  ex(1)+ex(3)/2 ex(1)       ex(1)];
y=[ex(2)+ex(4)/2 ex(2) ex(2)          ex(2)       ex(2)+ex(4)/2 ex(2)+ex(4)  ex(2)+ex(4)   ex(2)+ex(4) ex(2)+ex(4)/2];
case 2  % hexagon - bent sides
x=[ex(1)-ex(4)/6    ex(1)  ex(1)+ex(3)/2  ex(1)+ex(3) ex(1)+ex(3)+ex(4)/6   ex(1)+ex(3)  ex(1)+ex(3)/2   ex(1)         ex(1)-ex(4)/6];
y=[ex(2)+ex(4)/2    ex(2)  ex(2)          ex(2)       ex(2)+ex(4)/2         ex(2)+ex(4)  ex(2)+ex(4)     ex(2)+ex(4)   ex(2)+ex(4)/2];
case 3  % hexagon - flat sides
x=[ex(1)           ex(1)  ex(1)+ex(3)/2  ex(1)+ex(3)  ex(1)+ex(3)    ex(1)+ex(3)  ex(1)+ex(3)/2     ex(1)       ex(1)];
y=[ex(2)+ex(4)/2   ex(2)  ex(2)-ex(4)/5  ex(2)        ex(2)+ex(4)/2  ex(2)+ex(4)  ex(2)+ex(4)*6/5   ex(2)+ex(4) ex(2)+ex(4)/2];
case 4  % diamond
x=[ex(1)-ex(3)/2  ex(1)            ex(1)+ex(3)/2   ex(1)+ex(3)      ex(1)+1.5*ex(3)  ex(1)+ex(3)     ex(1)+ex(3)/2     ex(1)             ex(1)-ex(3)/2];
y=[ex(2)+ex(4)/2  ex(2)+ex(4)/12   ex(2)-ex(4)/3   ex(2)+ex(4)/12   ex(2)+ex(4)/2    ex(2)+7/8*ex(4) ex(2)+1.25*ex(4)  ex(2)+7/8*ex(4)   ex(2)+ex(4)/2];
case 5  % ellispe
x=linspace(-pi,pi,41);
y=sin(x)*ex(4)*.65+ex(2)+ex(4)/2;
x=cos(x)*ex(3)*.5+ex(1)+ex(3)/2;
case 6  % skewed hexagon
x=[ex(1)         ex(1)+ex(3)/7 ex(1)+ex(3)*4/7  ex(1)+ex(3) ex(1)+ex(3)   ex(1)+ex(3)*6/7  ex(1)+ex(3)*3/7 ex(1)       ex(1)];
y=[ex(2)+ex(4)/2 ex(2) ex(2)          ex(2)       ex(2)+ex(4)/2 ex(2)+ex(4)  ex(2)+ex(4)   ex(2)+ex(4) ex(2)+ex(4)/2];
case 7  % skewed hexagon
x=[ex(1)         ex(1) ex(1)+ex(3)*3/7  ex(1)+ex(3)*6/7 ex(1)+ex(3)   ex(1)+ex(3)  ex(1)+ex(3)*4/7 ex(1)+ex(3)/7       ex(1)];
y=[ex(2)+ex(4)/2 ex(2) ex(2)          ex(2)       ex(2)+ex(4)/2 ex(2)+ex(4)  ex(2)+ex(4)   ex(2)+ex(4) ex(2)+ex(4)/2];
end

function deletenode(h)
  cf=gcf;
  nodeinfo=get(h,'UserData');
  for i=1:length(nodeinfo.inarcs)
    deletearc(nodeinfo.inarcs(i))
  end
  for i=1:length(nodeinfo.outarcs)
    deletearc(nodeinfo.outarcs(i))
  end
  delete(nodeinfo.texthandle)
  figinfo=get(cf,'UserData');
  figinfo.nodes=figinfo.nodes(figinfo.nodes~=h);
  set(cf,'UserData',figinfo);
  delete(h)

function h=addarc(pp,pc,attachments)
  cf=gcf;
  figinfo=get(cf,'UserData');
  arcflag=figinfo.arcflag;
  loc1=min(1,getattachmentpoints(pp));
  loc2=min(1,getattachmentpoints(pc));
  if nargin>2
    loc1=loc1(attachments(1),:);
    loc2=loc2(attachments(2),:);
  else
    switch arcflag
      case 1
        attachments=[5 1]; loc1=loc1(5,:); loc2=loc2(1,:);
      case 2
        attachments=[3 7]; loc1=loc1(3,:); loc2=loc2(7,:);
      otherwise
        [loc1,loc2,attachments]=mindist(loc1,loc2);
    end
  end
  loc1=max(0,loc1);
  loc2=max(0,loc2);
  h=annotation('arrow',[loc1(1) loc2(1)],[loc1(2) loc2(2)]);
  set(h,'headwidth', get(h,'headwidth')*figinfo.archeadsize,...
        'headlength',get(h,'headlength')*figinfo.archeadsize,...
        'linewidth',figinfo.arcwidth)
  arcinfo=struct('tail',pp,'head',pc,'attachments',attachments);
  set(h,'UserData',arcinfo,'tag','arc')
  % add arc to parent and child node arclists
  nodeinfo=get(pp,'UserData');
  nodeinfo.outarcs=[nodeinfo.outarcs h];
  set(pp,'UserData',nodeinfo);
  nodeinfo=get(pc,'UserData');
  nodeinfo.inarcs=[nodeinfo.inarcs h];
  set(pc,'UserData',nodeinfo);
  set(h,'ButtonDownFcn',@mouseclick);
  uistack(get(h,'Parent'),'bottom')
  
function deletearc(h)
  cf=gcf;
  hh=findall(gca,'Tag','node');
  % remove arc from node arc lists
  for i=1:length(hh)
    nodeinfo=get(hh(i),'UserData');
    nodeinfo.inarcs=nodeinfo.inarcs(~ismember(nodeinfo.inarcs,h));
    nodeinfo.outarcs=nodeinfo.outarcs(~ismember(nodeinfo.outarcs,h));
    set(hh(i),'UserData',nodeinfo);
  end
  figinfo=get(cf,'UserData');
  figinfo.arcs=figinfo.arcs(figinfo.arcs~=h);
  set(cf,'UserData',figinfo);
  delete(h);
  

% find the attachment points for two graphic objects
% i.e., the best points to draw an arrow connecting the objects
function A=getattachmentpoints(h)
A=[get(h,'XData') get(h,'YData')];
switch size(A,1)
  case 41
    A=A(1:5:41,:); 
  case 9
    A=A(1:1:9,:); 
end

% mindist minimum distance between two sets of points
% USAGE
%   [loc1,loc2]=mindist(set1,set2);
% INPUTS
%   set1 : n1 x d matrix
%   set2 : n2 x d matrix
% OUTPUTS
%   loc1 : 1 x d point in set1 closest to set 2
%   loc2 : 1 x d point in set2 closest to set 1
function [loc1,loc2,attachments]=mindist(set1,set2)
n1=size(set1,1);
n2=size(set2,1);
mind=inf;
for i=1:n1
  for j=1:n2
    d=sum((set1(i,:)-set2(j,:)).^2);
    if d<mind, mind=d; attachments=[i j]; end
  end
end
loc1=set1(attachments(1),:);
loc2=set2(attachments(2),:);

% check if variables pc and pp are okay to join
% can't join if already joined or if they are the same
function ok=ok2addarc(pp,pc)
ok = true;
if pc==pp
  ok=false; return
end
nodecinfo=get(pc,'UserData');
nodepinfo=get(pp,'UserData');
if any(ismember(nodecinfo.inarcs,nodepinfo.outarcs))
  ok=false; return
end
if any(ismember(nodecinfo.outarcs,nodepinfo.inarcs))
  ok=false; return
end
if any(ismember(nodepinfo.inarcs,nodecinfo.outarcs))
  ok=false; return
end
if any(ismember(nodepinfo.outarcs,nodecinfo.inarcs))
  ok=false; return
end
figinfo=get(gcf,'UserData');
if figinfo.dag, ok=checkDAG(pp,pc); end
return

function DAG=checkDAG(hp,hc)
cf=gcf;
figinfo=get(cf,'UserData');
h=figinfo.nodes;
d=length(h);
A=zeros(d,d);
for i=1:d
  nodeinfo=get(h(i),'UserData');
  parents=nodeinfo.inarcs;
  for j=1:length(parents)
    pj=get(parents(j),'UserData');
    pj=pj.tail;
    A(pj==h,i)=1;
  end
end
A(hp==h,hc==h)=1;
[order,DAG]=topoorder(A); %#ok<ASGLU>
DAG=~DAG;

function redraw(h,event)
cf=gcf;
figinfo=get(cf,'UserData');
if isfield(figinfo,'nodes')
  for i=1:length(figinfo.nodes)
    hi=figinfo.nodes(i);
    nodeinfo=get(hi,'UserData');
    pos=get(nodeinfo.texthandle,'Position');
    pos(1)=ceil(pos(1)*512)/512;
    addnode(pos, get(nodeinfo.texthandle,'String'),...
          nodeinfo.type,nodeinfo.obs,nodeinfo.cpd,hi);
  end
  for i=1:length(figinfo.arcs)
    hi=figinfo.arcs(i);
    arcinfo=get(hi,'UserData');
    x=getattachmentpoints(arcinfo.tail);
    x=x(arcinfo.attachments(1),:);
    y=getattachmentpoints(arcinfo.head);
    y=y(arcinfo.attachments(2),:);
    set(hi,'X',[x(1) y(1)],'Y',[x(2) y(2)]);
  end
end

function shiftfig(dir)
cf=gcf;
figinfo=get(cf,'UserData');
if isfield(figinfo,'nodes')
  for i=1:length(figinfo.nodes)
    hi=figinfo.nodes(i);
    nodeinfo=get(hi,'UserData');
    pos=get(nodeinfo.texthandle,'Position');
    switch dir
      case 'up'    
        pos(2)=pos(2)+0.01;
      case 'down'
        pos(2)=pos(2)-0.01;
      case 'left'
        pos(1)=pos(1)-0.01;
      case 'right'
        pos(1)=pos(1)+0.01;
    end
    %pos(1)=ceil(pos(1)*512)/512;
    addnode(pos, get(nodeinfo.texthandle,'String'),...
          nodeinfo.type,nodeinfo.obs,nodeinfo.cpd,hi);
  end
  redraw;
end

function shrinkfig(dir)
cf=gcf;
figinfo=get(cf,'UserData');
if isfield(figinfo,'nodes')
  for i=1:length(figinfo.nodes)
    hi=figinfo.nodes(i);
    nodeinfo=get(hi,'UserData');
    pos=get(nodeinfo.texthandle,'Position');
    switch dir
      case 'shrink'
        pos(1)=pos(1)+0.05*(0.5-pos(1));
        pos(2)=pos(2)+0.05*(0.5-pos(2));
      case 'expand'
        pos(1)=pos(1)-0.05*(0.5-pos(1));
        pos(2)=pos(2)-0.05*(0.5-pos(2));
    end
    addnode(pos, get(nodeinfo.texthandle,'String'),...
          nodeinfo.type,nodeinfo.obs,nodeinfo.cpd,hi);
  end
  redraw;
end




%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks

function keypressfcn(src,evnt) %#ok<INUSL>
cf=gcf;
switch evnt.Key
  case 'delete'
    gui=get(gca,'UserData');
    figinfo=get(cf,'UserData');
    if ~isempty(gui)
      if gui.nodeselected
        yesnobox('Delete variable?',0.12,figinfo.fontname)   
        gui=get(gca,'UserData');          
        if gui.yesnoboxresult
          deletenode(gui.currenthandle);
          set(gca,'UserData',[]);
        end
      elseif gui.arcselected
        yesnobox('Delete arc?',[],figinfo.fontname)   
        gui=get(gca,'UserData');          
        if gui.yesnoboxresult
           deletearc(gui.currenthandle);
           set(gca,'UserData',[]);
        end
      end
    end
  case 'e'
    gui=get(gca,'UserData');
    if isempty(gui)
      %okbox('No variable selected')
    elseif gui.nodeselected
      %unselect(gui.currenthandle)
      figinfo=get(cf,'UserData');
      figinfo.nodeedit(gui.currenthandle,@closenodebox)
      set(gca,'UserData',[]);
      redraw
    end
  case 'f1'
    helpbox
  case 'f2'
    figinfo=get(cf,'UserData');
    color=uisetcolor('Choose node color');
    if length(color)==3
      figinfo.nodecolor=color;
      set(cf,'UserData',figinfo);
      redraw
    end
  case 'f3'
    figinfo=get(cf,'UserData');
    color=uisetcolor('Choose background color');
    if length(color)==3
      figinfo.backgroundcolor=color;
      set(cf,'Color',color)
      set(cf,'UserData',figinfo);
    end
  case 'f4'
    figinfo=get(cf,'UserData');
    fontname=uisetfont('Select font name');
    if ~isnumeric(fontname) && ~strcmp(fontname.FontName,figinfo.fontname)
      figinfo.fontname=fontname.FontName;
      set(cf,'UserData',figinfo);
      redraw
    end
  case 'f5'
    if strcmp(get(cf,'Toolbar'),'none')
      set(cf,'ToolBar','figure','MenuBar','figure')
    else
      set(cf,'ToolBar','none','MenuBar','none')
    end
  case 'f6'
    dname = inputdlg('Enter diagram name:','',1,{'D'},'on');
    if ~isempty(dname)
      figinfo=get(cf,'UserData');
      typestring=figinfo.typestring;
      [~,ii]=unique(figinfo.types); 
      updatediagram([],gcf,2,dname{1},typestring(ii));
      %updatediagram([],gcf,2,dname{1},typestring);
    end
  case 'f7'
    figinfo=get(cf,'UserData');
    if figinfo.confirmationboxes
      figinfo.confirmationboxes=false;
    else
      figinfo.confirmationboxes=true;
    end
    set(cf,'UserData',figinfo);
  case 'f8'
    figinfo=get(cf,'UserData');
    figinfo.boxsize=figinfo.boxsize/1.05;
    set(cf,'UserData',figinfo);
  case 'f9'
    figinfo=get(cf,'UserData');
    figinfo.boxsize=figinfo.boxsize*1.05;
    set(cf,'UserData',figinfo);
  case 'pagedown'  
    figinfo=get(cf,'UserData');
    figinfo.fontsize=figinfo.fontsize/1.1;
    set(cf,'UserData',figinfo);
    redraw
  case 'pageup'  
    figinfo=get(cf,'UserData');
    figinfo.fontsize=figinfo.fontsize*1.1;
    set(cf,'UserData',figinfo);
    redraw
  case 'uparrow'
    shiftfig('up')
  case 'downarrow'
    shiftfig('down')
  case 'leftarrow'
    shiftfig('left')
  case 'rightarrow'
    shiftfig('right')    
  case 'home'
    shrinkfig('expand')    
  case 'end'
    shrinkfig('shrink')
end

function select(h)
cf=gcf;
figinfo=get(cf,'UserData');
switch get(h,'Tag')
  case 'node'
    nodedata=get(h,'UserData');
    if isfield(nodedata,'cpd') && isnumeric(nodedata.cpd) && numel(nodedata.cpd)==1 && nodedata.cpd==2
      set(h,'EdgeColor',[1 0 0],'LineWidth',figinfo.nodeedgewidth*3)
    else
      set(h,'EdgeColor',[1 0 0],'LineWidth',figinfo.nodeedgewidth*3)
    end
  case 'arc'
      set(h,'Color',[1 0 0],'LineWidth',figinfo.arcwidth*2)
end

function unselect(h)
cf=gcf;
figinfo=get(cf,'UserData');
switch get(h,'Tag')
  case 'node'
    nodedata=get(h,'UserData');
    if isnumeric(nodedata.cpd) && ~isempty(nodedata.cpd) && nodedata.cpd==2
      set(h,'EdgeColor',[0 0 0],'LineWidth',figinfo.nodeedgewidth)
    else
      set(h,'EdgeColor',[0 0 0],'LineWidth',figinfo.nodeedgewidth)
    end
  case 'arc'
      set(h,'Color',[0 0 0],'LineWidth',figinfo.arcwidth)
end
set(gca,'UserData',[]);

function mouseclick(h,event) 
cf=gcf;
ca=gca;
hf=findall(0,'Tag','NodeDefinitionWindow');
if ~isempty(hf), close(hf); end
hf=findall(0,'Tag','YesNoBox');
if ~isempty(hf), close(hf); end
hf=findall(0,'Tag','OKBox');
if ~isempty(hf), close(hf); end
type=get(h,'Tag');
% text clicked - change handle to patch
if strcmp(type,'nodetext')
  h=get(h,'UserData');
  type='node';
end
seltype=get(cf,'SelectionType');
gui=get(ca,'UserData');
figinfo=get(cf,'UserData');

% turn off selection
if isstruct(gui)&& gui.currenthandle==h
  switch type
    case 'node'
      if all(get(gui.currenthandle,'EdgeColor')==[1 0 0]) 
        unselect(gui.currenthandle)
        set(ca,'UserData',[]);
        if strcmp(seltype,'alt'), return; end
      end
    case 'arc'
      if all(get(gui.currenthandle,'Color')==[1 0 0])
        unselect(gui.currenthandle)
        set(ca,'UserData',[]);
        if strcmp(seltype,'alt'), return; end
      end
  end
  gui=[];
end

switch type
  case 'node'
    switch seltype
      case 'alt'
        if isstruct(gui)
          if gui.nodeselected
            if ok2addarc(gui.currenthandle,h)
              nodedata=get(gui.currenthandle,'UserData');
              qstring={'Add an arc from',...
                 ['\it ' get(nodedata.texthandle,'String')],'to',[]};
              nodedata=get(h,'UserData');
              qstring{4}=['\it ' get(nodedata.texthandle,'String') '?'];
              yesnobox(qstring,[],figinfo.fontname)   
              gui=get(ca,'UserData'); 
              if ~isfield(gui,'yesnoboxresult') 
                if isfield(gui,'currenthandle')
                  unselect(gui.currenthandle)
                end
                set(ca,'UserData',[])
                set(cf,'UserData',figinfo);
              elseif gui.yesnoboxresult
                figinfo.arcs(end+1)=addarc(gui.currenthandle,h);
                if isfield(gui,'currenthandle')
                  unselect(gui.currenthandle)
                end
                set(ca,'UserData',[])
                set(cf,'UserData',figinfo);
              end
            else
              okbox('Invalid arc selection',[],figinfo.fontname)
            end
          else
            unselect(gui.currenthandle)
            select(h)
            gui.currenthandle=h;
            gui.nodeselected=true;
            gui.arcselected=false;
            set(ca,'UserData',gui);
          end
        else
          select(h)
          gui.currenthandle=h;
          gui.nodeselected=true;
          gui.arcselected=false;
          set(ca,'UserData',gui)
        end 
      otherwise
        if isstruct(gui)
          unselect(gui.currenthandle)
          set(ca,'UserData',[]);
        end
        startmovenode(h);
    end
  case 'arc'
    switch seltype
      case 'alt' 
        select(h)
        gui=get(ca,'UserData');
        if isstruct(gui)
          unselect(gui.currenthandle);
          gui=[]; 
        end
        gui.currenthandle=h;
        gui.arcselected=true;
        gui.nodeselected=false;
        set(ca,'UserData',gui)
      otherwise
        if isstruct(gui)
          unselect(gui.currenthandle)
          set(ca,'UserData',[]);
        end
        startmovearc(h);
    end
  otherwise
    switch seltype
      case 'alt'
        figinfo=get(cf,'UserData');
        figinfo.nodeedit([],@closenodebox)
      case 'normal'
    end
end


function startmovenode(src)
  cf=gcf;
% Unpack gui object
  gui = get(gca,'UserData');
  if ~isempty(gui)
    stopmovenode(gui.currenthandle)
  end
  % Remove mouse pointer
  %set(gcf,'PointerShapeCData',nan(16,16));
  %set(gcf,'Pointer','custom');
  % Set callbacks
  gui.currenthandle=src;
  % Store starting point of the object
  gui.startpoint = get(cf,'CurrentPoint');
  nodedata=get(src,'UserData');
  % find the associated text box
  uistack(src,'top')
  uistack(nodedata.texthandle,'top')
  %uistack(nodedata.texthandle,'bottom')
  %uistack(src,'bottom')
  gui.texthandle=nodedata.texthandle;
  %gui.corners=getattachmentpoints(src);
  gui.textpos=get(gui.texthandle,'Position');
  gui.XData=get(src,'Xdata');
  gui.YData=get(src,'Ydata');
  gui.inarcs=nodedata.inarcs;
  gui.inlocs=zeros(length(nodedata.inarcs),4);
  for i=1:length(nodedata.inarcs)
    gui.inlocs(i,:)=[get(gui.inarcs(i),'X') get(gui.inarcs(i),'Y')];
  end
  gui.outarcs=nodedata.outarcs;
  gui.outlocs=zeros(length(nodedata.outarcs),4);
  for i=1:length(gui.outarcs)
    gui.outlocs(i,:)=[get(gui.outarcs(i),'X') get(gui.outarcs(i),'Y')];
  end
  thisfig = gcbf();
  set(thisfig,'WindowButtonMotionFcn',@movenode);
  set(thisfig,'WindowButtonUpFcn',@stopmovenode);
  % Store gui object
  set(gca,'UserData',gui);




function movenode(src,evnt) 
cf=gcf;
% Unpack gui object
gui = get(gca,'UserData');
if ~isfield(gui,'startpoint') || isempty(gui.startpoint)
  return
end
% Do "smart" positioning of the object, relative to starting point...
try
  pos = get(cf,'CurrentPoint')-gui.startpoint;
  XYData{1}=gui.XData + pos(1,1);
  XYData{2}=gui.YData + pos(1,2);
  % only move if object stays inside the figure
  if all(XYData{1}<=1) && all(XYData{2}<=1) && all(XYData{1}>=0) && all(XYData{2}>=0)
    for i=1:length(gui.inarcs)
      XY=gui.inlocs(i,:);
      set(gui.inarcs(i),'X',XY(1:2)+[0 pos(1)],'Y',XY(3:4)+[0 pos(2)])
    end
    for i=1:length(gui.outarcs)
      XY=gui.outlocs(i,:);
      set(gui.outarcs(i),'X',XY(1:2)+[pos(1) 0],'Y',XY(3:4)+[pos(2) 0])
    end
    set(gui.currenthandle,'XData',XYData{1});
    set(gui.currenthandle,'YData',XYData{2});
    set(gui.texthandle,'Position',gui.textpos+[pos 0])
    drawnow;
  end
catch
  stopmovenode
  drawnow;
end


function stopmovenode(src,evnt) 
% Clean up the evidence ...
thisfig = gcbf();
set(thisfig,'Pointer','arrow');
set(thisfig,'WindowButtonUpFcn','');
set(thisfig,'WindowButtonMotionFcn','');
gui = get(gca,'UserData');
if isfield(gui,'currenthandle')
  unselect(gui.currenthandle)
end
set(gca,'UserData',[]);



function startmovearc(h,event) %#ok<*INUSD>
% Remove mouse pointer
%set(gcf,'PointerShapeCData',nan(16,16));
%set(gcf,'Pointer','custom');

gui = get(gca,'UserData');
if ~isempty(gui)
  stopmovearc(gui.currenthandle)
end

% get starting point of the object
cf=gcf;
startpoint = get(cf,'CurrentPoint');
X=get(h,'X');
Y=get(h,'Y');
arcinfo=get(h,'UserData');
% determine if head or tail moves
if abs(X(1)-startpoint(1))+abs(Y(1)-startpoint(2)) < ...
   abs(X(2)-startpoint(1))+abs(Y(2)-startpoint(2))
   gui.headmoves=false;
   gui.movenode=arcinfo.tail;
else
  gui.headmoves=true;
  gui.movenode=arcinfo.head;
end
gui.currenthandle=h;
gui.startpoint = startpoint; 
gui.arcselected=true;
gui.X=X;
gui.Y=Y;
% Store gui object
set(gca,'UserData',gui);

thisfig = gcbf();
set(thisfig,'WindowButtonMotionFcn',@movearc);
set(thisfig,'WindowButtonUpFcn',@stopmovearc);


function movearc(src,evnt)
% Unpack gui object
gui = get(gca,'UserData');
if ~isfield(gui,'startpoint') || isempty(gui.startpoint)
  return
end
% Do "smart" positioning of the object, relative to starting point...
pos = get(gcf,'CurrentPoint')-gui.startpoint;
if gui.headmoves
  set(gui.currenthandle,'X',gui.X + [0 pos(1,1)]);
  set(gui.currenthandle,'Y',gui.Y + [0 pos(1,2)]);
else
  set(gui.currenthandle,'X',gui.X + [pos(1,1) 0]);
  set(gui.currenthandle,'Y',gui.Y + [pos(1,2) 0]);
end
drawnow;


function stopmovearc(src,evnt)
% Clean up the evidence ...
thisfig = gcbf();
gui = get(gca,'UserData');
set(gcf,'Pointer','arrow');
set(thisfig,'WindowButtonUpFcn','');
set(thisfig,'WindowButtonMotionFcn','');
arcinfo=get(gui.currenthandle,'UserData');
A=getattachmentpoints(gui.movenode);
pos=[get(gui.currenthandle,'X') get(gui.currenthandle,'Y')];
if gui.headmoves
  [tmp,pos1]=min(sum(abs(bsxfun(@minus,A,pos([2 4]))),2));  
  arcinfo.attachments(2)=pos1;
  pos1=A(pos1,:);
  pos2=gui.X;
  set(gui.currenthandle,'X',[pos2(1) pos1(1)]);
  pos2=gui.Y;
  set(gui.currenthandle,'Y',[pos2(1) pos1(2)]);
else
  [tmp,pos1]=min(sum(abs(bsxfun(@minus,A,pos([1 3]))),2)); 
  arcinfo.attachments(1)=pos1;
  pos1=A(pos1,:);
  pos2=gui.X;
  set(gui.currenthandle,'X',[pos1(1) pos2(2)]);
  pos2=gui.Y;
  set(gui.currenthandle,'Y',[pos1(2) pos2(2)]);
end
set(gui.currenthandle,'UserData',arcinfo);
set(gca,'UserData',[]);
drawnow;

function helpbox
Message={...
  'Left click on a blank area to add a variable there',  ...
  'Right click to select or deselect a variable or arc',  ...
  'Right click a second variable to add an arc from the selected variable', ...
  'Delete key to remove a selected variable or arc',  ...
  'e key to edit a selected variable',...
  'Left click and hold to move a variable or arc',  ...
  'Arcs can only be moved to one of 8 attachment points on a variable; ',  ...
  '     a moved arc will snap to the nearest attachment point',  ...
  '  ', ...
  'F2 to change variable color',  ...
  'F3 to change background color',  ...
  'F4 to change the font',  ...
  'F5 toggle the Toolbar and Menubar on and off',...
  'F6 to print out current diagram information (runs updatediagram)',...
  'F7 to toggle confirmation boxes',...
  'Arrow, PgUp, PgDn, Home and End keys alter size and location'};
okbox(Message,[.3 .45],[],'Left')


%%%%%%%%%% UI control for new variables
function getnewvariable(h,closebox)
cf=gcf;
figinfo=get(cf,'UserData');
if isempty(h)
  name='';
  type=1;
  obs=1;
  cpd='';
else
  nodeinfo=get(h,'UserData');
  name=get(nodeinfo.texthandle,'String');
  type=nodeinfo.type;
  obs=nodeinfo.obs;
  cpd=nodeinfo.cpd;
end
pos=get(cf,'CurrentPoint');
data={pos,h};
q=get(cf,'Position');
pos(1)=pos(1)*q(3)+q(1);
pos(2)=pos(2)*q(4)+q(2);
pos=[pos 0.3*figinfo.boxsize 0.4*figinfo.boxsize];
% open towards the center from the current point
if data{1}(1)>0.5, pos(1)=pos(1)-pos(3); end
if data{1}(2)>0.5, pos(2)=pos(2)-pos(4); end
% Create and then hide the GUI as it is being constructed.
f = figure('Visible','off');
set(f,'Units','normalized','OuterPosition',pos,...
    'NumberTitle','off',...
    'Toolbar','none','MenuBar','none','Color',[0.85 0.85 0.85],...
    'Tag','NodeDefinitionWindow',...
    'Interruptible','off','DefaultUIControlFontSize',10);  
  
backgroundcolor=get(f,'Color');
  
if isempty(h)
  set(f,'Name','Define new variable')
else
  set(f,'Name','Edit existing variable')
end
  
uicontrol('Parent',f,'Style','text','Units','normalized',...
  'Position',[.05 .845 .3 .1], 'String','Variable name:  ','Visible','on','HorizontalAlignment','right',...
   'BackgroundColor',[0.85 0.85 0.85]);
uname = uicontrol('Parent',f,'Style','edit','Units','normalized',...
  'Position',[.35 .865 .6 .1], 'String',name,'Visible','on','HorizontalAlignment','left',...
  'Callback',@changename);

hbt=uibuttongroup('Parent',f,'FontSize',10,'Title','Variable type',...
  'Position',[.05 .2 .4 .625],'BackgroundColor',backgroundcolor,'SelectionChangeFcn',@changestyle);
ht=zeros(6,1);
ht(1)=uicontrol(hbt,'Style','radiobutton','Units','normalized',...
  'String','state (s)',       'Value',0,'Position',[.05 11/13 .8 .12],'BackgroundColor',backgroundcolor);
ht(2)=uicontrol(hbt,'Style','radiobutton','Units','normalized',...
  'String','action (a,d)',    'Value',0,'Position',[.05 9/13 .8 .12],'BackgroundColor',backgroundcolor);
ht(3)=uicontrol(hbt,'Style','radiobutton','Units','normalized',...
  'String','future state (f)','Value',0,'Position',[.05 7/13 .8 .12],'BackgroundColor',backgroundcolor);
ht(4)=uicontrol(hbt,'Style','radiobutton','Units','normalized',...
  'String','utility (u,r)',   'Value',0,'Position',[.05 5/13 .8 .12],'BackgroundColor',backgroundcolor);
ht(5)=uicontrol(hbt,'Style','radiobutton','Units','normalized',...
  'String','chance (c)',      'Value',0,'Position',[.05 3/13 .8 .12],'BackgroundColor',backgroundcolor);
ht(6)=uicontrol(hbt,'Style','radiobutton','Units','normalized',...
  'String','parameter (p)',   'Value',0,'Position',[.05 1/13 .8 .12],'BackgroundColor',backgroundcolor);
set(hbt,'SelectedObject',ht(type));
 

hbo=uibuttongroup('Parent',f,'FontSize',10,'Title','Observed',...
  'Position',[.55 .575 .4 .25],'BackgroundColor',backgroundcolor,'SelectionChangeFcn',@changeobs);
ho=zeros(2,1);
ho(1)=uicontrol(hbo,'Style','radiobutton','Units','normalized',...
  'String','yes',       'Value',0,'Position',[.05 3/5 .8 .4],'BackgroundColor',backgroundcolor);
ho(2)=uicontrol(hbo,'Style','radiobutton','Units','normalized',...
  'String','no',       'Value',0,'Position',[.05 1/5 .8 .4],'BackgroundColor',backgroundcolor);
if obs
  set(hbo,'SelectedObject',ho(1));
else
  set(hbo,'SelectedObject',ho(2));
end   

uicontrol('Parent',f,'Style','text','Units','normalized',...
  'Position',[.55 .35 .4 .15], 'String','Conditional Probability Distribution (CPD):  ','Visible','on','HorizontalAlignment','left',...
   'BackgroundColor',[0.85 0.85 0.85]);
ucpd = uicontrol('Parent',f,'Style','edit','Units','normalized',...
  'Position',[.55 .25 .4 .1], 'String',cpd,'Visible','on','HorizontalAlignment','left',...
  'Callback',@changecpd);


uicontrol(f,'Style','pushbutton','Units','normalized',...
  'String','Accept',       'Value',0,'Position',[.25 .025 .2 .125],...
  'Callback',@closenodebox,'UserData','accept');
uicontrol(f,'Style','pushbutton','Units','normalized',...
  'String','Cancel',       'Value',0,'Position',[.55 .025 .2 .125],...
  'Callback',@closenodebox,'UserData','cancel');

data=[{name,type,obs,cpd} data];



set(f,'UserData',data);
set(gca,'Visible','off')
%Make the GUI visible.
drawnow
set(f,'Visible','on')
drawnow
uicontrol(uname)

function changename(src,evnt)
cf=gcf;
data=get(cf,'UserData');
s=get(src,'String');
if isempty(s)
  uicontrol(src)
else
  data{1}=s;
  set(cf,'UserData',data);
end



function changestyle(src,evnt)
cf=gcf;
data=get(cf,'UserData');
% variable type
data{2}=7-find(ismember(get(src,'Children'), get(src,'SelectedObject')));
set(cf,'UserData',data);

function changeobs(src,evnt)
cf=gcf;
data=get(cf,'UserData');
ii=find(ismember(get(src,'Children'), get(src,'SelectedObject')));
switch ii
  case 1
    data{3}=0;
  case 2
    data{3}=1;
end
set(cf,'UserData',data);

function changecpd(src,evnt)
cf=gcf;
data=get(cf,'UserData');
s=get(src,'String');
data{4}=s;
set(cf,'UserData',data);

function closenodebox(src,evnt)
cf=gcf;
data=get(cf,'UserData');
action=get(src,'UserData');
switch action
  case 'accept'
    if isempty(data{1})
      msgbox('Must enter a valid variable name', 'Warning','warn','modal');
    else
      close(cf)
      if isempty(data{end})
        figinfo=get(gcf,'UserData');
        figinfo.nodes(end+1)=addnode(data{5},data{1:4});
        set(gcf,'UserData',figinfo)
      else
        nodeinfo=get(data{end},'UserData');
        loc=get(nodeinfo.texthandle,'Position');
        addnode(loc,data{1:4},data{end});
        unselect(data{end})
      end
    end
  case 'cancel'
    close(cf)
end


function yesnobox(Question,size,fontname)
cf=gcf;
figinfo=get(cf,'UserData');
gui=get(gca,'UserData');
if ~figinfo.confirmationboxes
  gui.yesnoboxresult=true;
  set(gca,'UserData',gui);
  return
end
if ~iscell(Question)
  Question={Question};
end

nblocks=length(Question)+3;
if nargin<2
  size=[0.1 0.03*nblocks];
end

factor=0.05;
if nargin<2 || isempty(size)
  size=[0.1 factor*nblocks];
elseif length(size)==1
  size=[size factor*nblocks];
end

if nargin<3 || isempty(fontname)
  fontname='Times New Roman';
end

gui.yesnoboxresult=false;
set(gca,'UserData',gui);

units=get(gcf,'Units');
set(cf,'Units','Normalized');
pos=get(cf,'CurrentPoint');
q=get(cf,'Position');
set(cf,'Units',units);
figpos=pos;
pos(1)=pos(1)*q(3)+q(1);
pos(2)=pos(2)*q(4)+q(2);
pos=[pos size*figinfo.boxsize];
% open towards the center from the current point
if figpos(1)>0.5, pos(1)=pos(1)-pos(3); end
if figpos(2)>0.5, pos(2)=pos(2)-pos(4); end
color=[0.85 0.85 0.85];
% Create and then hide the GUI as it is being constructed.
f = figure('Visible','off');
set(f,'Units','normalized','OuterPosition',pos,...
    'NumberTitle','off','Name','',...
    'Toolbar','none','MenuBar','none','Color',color,...
    'Tag','YesNoBox','Interruptible','off');  

for i=1:length(Question)  
  if length(Question{i})>4 && strcmp(Question{i}(1:4),'\it ')
    FontAngle='italic';
    Question{i}=Question{i}(5:end);
  else
    FontAngle='normal';
  end
  uicontrol('Parent',f,'Style','text','Units','normalized',...
    'Position',[.05 1-(i+1)/nblocks .9 1/nblocks], 'String',Question{i},'Visible','on',...
    'HorizontalAlignment','center',...
    'BackgroundColor',color,'FontSize',12,'FontAngle',FontAngle,...
    'FontName',fontname);
end

yesbutton=uicontrol(f,'Style','pushbutton','Units','normalized',...
  'String','Yes',       'Value',0,'Position',[.1 0.5/nblocks .3 1/nblocks],...
  'Callback',@closeyesnobox,'UserData','Yes');
uicontrol(f,'Style','pushbutton','Units','normalized',...
  'String','No',       'Value',0,'Position',[.6 0.5/nblocks .3 1/nblocks],...
  'Callback',@closeyesnobox,'UserData','No');

%Make the GUI visible.
drawnow
set(f,'Visible','on')
uicontrol(yesbutton)
uiwait(f)

function closeyesnobox(src,evnt)
action=get(src,'UserData');
close(gcf)
gui=get(gca,'UserData');
switch action
  case 'Yes'
    gui.yesnoboxresult=true;
  otherwise
    gui.yesnoboxresult=false;
end
set(gca,'UserData',gui);


function okbox(Message,size,fontname,HorizontalAlignment)
cf=gcf;
figinfo=get(cf,'UserData');
gui=get(gca,'UserData');
if ~iscell(Message)
  Message={Message};
end
nblocks=length(Message)+3;
factor=0.05;
if nargin<2 || isempty(size)
  size=[0.1 factor*nblocks];
elseif length(size)==1
  size=[size factor*nblocks];
end
if nargin<3 || isempty(fontname)
  fontname='Times New Roman';
end
if nargin<4 || isempty(HorizontalAlignment)
  HorizontalAlignment='center';
end
units=get(cf,'Units');
set(gcf,'Units','Normalized');
pos=get(cf,'CurrentPoint');
q=get(cf,'Position');
set(cf,'Units',units);
figpos=pos;
pos(1)=pos(1)*q(3)+q(1);
pos(2)=pos(2)*q(4)+q(2);
pos=[pos size*figinfo.boxsize];
% open towards the center from the current point
if figpos(1)>0.5, pos(1)=pos(1)-pos(3); end
if figpos(2)>0.5, pos(2)=pos(2)-pos(4); end
color=[0.85 0.85 0.85];
% Create and then hide the GUI as it is being constructed.
f = figure('Visible','off');
set(f,'Units','normalized','OuterPosition',pos,...
    'NumberTitle','off','Name','',...
    'Toolbar','none','MenuBar','none','Color',color,...
    'Tag','OKBox','Interruptible','off');  
for i=1:length(Message)  
  if length(Message{i})>4 && strcmp(Message{i}(1:4),'\it ')
    FontAngle='italic';
    Message{i}=Message{i}(5:end);
  else
    FontAngle='normal';
  end  
  uicontrol('Parent',f,'Style','text','Units','normalized',...
    'Position',[.05 1-(i+1)/nblocks .9 1/nblocks], 'String',Message{i},'Visible','on',...
    'BackgroundColor',color,'FontSize',12,'FontAngle',FontAngle,...
    'HorizontalAlignment',HorizontalAlignment,...
    'FontName',fontname);
end
okbutton=uicontrol(f,'Style','pushbutton','Units','normalized',...
  'String','OK',       'Value',0,'Position',[.4 0.5/nblocks .3 1/nblocks],...
  'Callback',@(src,evnt)close(gcf),'UserData','OK');
%Make the GUI visible.
drawnow
set(f,'Visible','on')
uicontrol(okbutton)
uiwait(f)


% topoorder Topological ordering for a Directed Acyclic Graph (DAG)
% USAGE
%   [L,errflag]=topoorder(A);
% INPUT
%   A : an adjacency matrix or an influence diagram structure
% OUTPUTS
%   L : a vector with a variable ordering
%   errflag : 0-successful 1-graph not a DAG (returns empty L)
function [L,errflag]=topoorder(A)
if isstruct(A)
  A=adjacency(A);
end
errflag=0;
n=size(A,1);
k=0;
L=zeros(1,n);
S=find(sum(A,1)==0);
while 1
  if isempty(S), break; end
  next=S(1); S(1)=[];
  k=k+1; L(k)=next;
  children=find(A(next,:)~=0);
  A(next,:)=0;
  for i=1:length(children)
    if all(A(:,children(i))==0)
      S=[S children(i)];
    end
  end
end
if any(any(A))
  L=[];
  if nargout<2
    warning('graph is not acyclic'); 
  else
    errflag=1;
  end
end

