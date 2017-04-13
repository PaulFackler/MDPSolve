% mdpplot Plots results from MDPSOLVE
% USAGE
%   h=mdpplot(S,V,ind,labels,options);
% INPUTS
%   S : nxd matrix of state variable values
%   V : nx1 vector of values to be plotted
%   ind : vector of indices of variables to be plotted
%            ind(1) : x-axis
%            ind(2) : y-axis
%            ind(3) : subplot rows (one row if ind(3)=0)
%            ind(4) : subplot columns (one column if ind(4)=0)
%   labels : 4-element cell array with labels for the 4 variables
%              (set to empty if no labels are desired)
%   options : a structure variable with the following allowable fields
%               edges         : 1 to plot edge borders around each cell
%               squareplot    : 1 to make the x and y axes of equal size
%               addlegend     : 1 to add a legend 
%               vertical      : 1 to make the legend have vertical orientation
%                                 (placed outside to the east - if 0 placed to south)
%               legendlabels  : a cell array of labels for a legend
%                                 there must be exactly as many as there are 
%                                 unique values of V 
%               figuretitle   : title placed in the figure header
%               legendtitle   : string to put over the legend box
%               colorbartype  : 1 to make the legend be a colorbar type rather than discrete
%               grayscale     : 1 to use gray color scheme
%               noticklabels  : 1 to supress ticklabels on x and y axes
% OUTPUTS
%   h  : a matrix of plot handles, one for each subplot
%   hl : handle for the legend

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

function [h,hl,hc]=mdpplot(S,V,ind,labels,options)  
  hl=[];
  if nargin<5, options=[]; end
  if nargin<4, labels={};  end
  getopts(options, ...
      'clim',         [],...
      'edges',        0, ...     % 1 if each cell is framed
      'squareplot',   0, ...
      'addlegend',    0, ...
      'vertical',     0, ...     % 1 if legend is vertical
      'legendlabels', {}, ...
      'colorbartype', 0, ...
      'grayscale',    0, ...
      'facecolortype', 'flat', ...
      'figuretitle',   '', ...
      'legendtitle',   '', ...
      'fontsize',      get(gca,'fontsize'), ...
      'fontname',      'Times New Roman', ...
      'noticklabels', 0);  
  
  set(gcf,'Name',figuretitle,'color',[1 1 1]);
    
  if isnumeric(ind)
    ii=ind;
    ind=cell(1,length(ii));
    for i=1:length(ind), ind{i}=ii(i); end
    clear ii
  end
 
  % set edge color for patches
  if edges
    edgecolor=[0 0 0];
  else
    edgecolor='none';
  end
  
  % define a gray color map
  % avoids very black and very white
  if grayscale
    colormap('default')
    colormap('gray'); 
    C=colormap; 
    colormap(flipud(C(9:56,:)));
  end
  
  if addlegend 
    if vertical
      hlegheight=0;
      vlegwidth=0.1;
    else
      hlegheight=0.1;
      vlegwidth=0;
    end
  else
    hlegheight=0;
    vlegwidth=0;
  end
  
  xlabheight = 0.025;
  ylabwidth  = 0.025;
  
  
  %xlabheight = 0;
  %ylabwidth  = 0;
   
  if length(ind)>=3 && all(ind{3}>0) 
    zlabels=true;
    zlabheight = 0.05;
  else
    zlabels=false;
    zlabheight = 0;
  end
  
   if length(ind)>=4 && all(ind{4}>0) 
    wlabels=true;
    wlabwidth = 0.05;
  else
    wlabels=false;
    wlabwidth = 0;
  end
  
  
  % gets the unique elements of S1 and S2
  % defines patch locations
  [s1,temp,xind]=unique(S(:,ind{1}));
  d1=(s1(2)-s1(1))/2;
  [s2,temp,yind]=unique(S(:,ind{2}));
  d2=(s2(2)-s2(1))/2;
  X=[s1'-d1;s1'+d1;s1'+d1;s1'-d1];  % x-axis patch coordinates
  Y=[s2'-d2;s2'-d2;s2'+d2;s2'+d2];  % y-axis patch coordinates
  abx=[s1(1)-d1 s1(end)+d1];        % x-axis min and max values
  aby=[s2(1)-d2 s2(end)+d2];        % y-axis min and max values
  
  % get Z-axis information
  vvals=unique(V);
  vmin=vvals(1);
  vmax=vvals(end);
  Z=ones(4,1)*V(:)';                % z-axis patch coordinates
  if isempty(clim)
  if vmax>vmin
    clim=[vmin,vmax]; 
  elseif addlegend && (exist('leglabels','var') && length(leglabels)>1)
    error('Cannot determine scaling for legend')
  else
    clim=[vmin,vmin+1];
  end
  end
  
  % gets unique values and number of values for 4th variable
  if length(ind)>=4 
    [Sr,temp,ri]=unique(S(:,ind{4}),'rows');
    nrows=size(Sr,1);
  else
    ri=ones(size(S,1),1);
    nrows=1;
  end
  % gets unique values and number of values for 3th variable
  if length(ind)>=3 && all(ind{3}>0) 
    Sc=unique(S(:,ind{3}),'rows');
    ncols=size(Sc,1);
  else
    ncols=1;
  end
  
  clf
  hfigure=gcf;
  hpanels=zeros(7,1);
  % graph panel - put first in order so it does not cover other panels
  hpanels(7)= uipanel('Parent',hfigure, 'BackgroundColor',[1 1 1], 'borderwidth', 0, ...
                    'Position',[ylabwidth+wlabwidth hlegheight+xlabheight 1-vlegwidth-ylabwidth-wlabwidth 1-hlegheight-xlabheight-zlabheight]);
  % horizontal legend panel placed at bottom of the figure
  if hlegheight>0
    %hpanels(1) = uipanel('Parent',hfigure, 'BackgroundColor',[1 1 1], 'borderwidth', 0, ...
    %                   'Position',[ylabwidth+wlabwidth 0 1-ylabwidth-wlabwidth hlegheight]);
    hpanels(1) = uipanel('Parent',hfigure, 'BackgroundColor',[1 1 1], 'borderwidth', 0, ...
                         'Position',[ylabwidth 0 1-ylabwidth hlegheight]);
  end
  % vertical legend panel placed at right side of the figure
  if vlegwidth>0
    hpanels(2) = uipanel('Parent',hfigure, 'BackgroundColor',[1 1 1], 'borderwidth', 0, ...
                       'Position',[1-vlegwidth xlabheight vlegwidth 1-xlabheight-zlabheight]);
  end
  if nrows>1 && ncols>1
  % panels for x and y labels
    hpanels(3)= uipanel('Parent',hfigure, 'BackgroundColor',[1 1 1], 'borderwidth', 0, ...
                        'Position',[ylabwidth+wlabwidth hlegheight 1-vlegwidth-ylabwidth-wlabwidth xlabheight]);
    hpanels(4)= uipanel('Parent',hfigure, 'BackgroundColor',[1 1 1], 'borderwidth', 0, ...
                        'Position',[wlabwidth hlegheight+xlabheight ylabwidth 1-hlegheight-xlabheight-zlabheight]);
  end
  % panel for z labels at the top of the plot
  if zlabels
    hpanels(5)= uipanel('Parent',hfigure, 'BackgroundColor',[1 1 1], 'borderwidth', 0, ...
                      'Position',[ylabwidth+wlabwidth 1-zlabheight 1-vlegwidth-ylabwidth-wlabwidth zlabheight]);
  end
  % panel for the w labels at the left side of the plot
  if wlabels
    hpanels(6)= uipanel('Parent',hfigure, 'BackgroundColor',[1 1 1], 'borderwidth', 0, ...
                      'Position',[0 hlegheight+xlabheight wlabwidth 1-hlegheight-xlabheight-zlabheight]);
  end
                  
 gwidth=1/ncols;
 gheight=1/nrows;

  
  h=zeros(nrows,ncols); % handles to patch objects
  for i=1:nrows
    for j=1:ncols
      h(i,j)=axes('parent',hpanels(7),'outerposition',[(j-1)*gwidth (i-1)*gheight gwidth gheight],'Color',[1 1 1]);
    end
  end
  
  % loop over variables 3 (j) and 4 (i)
  for i=1:nrows
    indi=find(ri==i);
    Si=S(indi,:);
    if length(ind)>=3 && all(ind{3}>0) 
      [Sc,temp,cj]=unique(Si(:,ind{3}),'rows');
    else
      cj=ones(length(indi),1);
    end
    for j=1:ncols
      indij=(indi(cj==j));
      %axes(h(i,j))
      set(gcf,'CurrentAxes',h(i,j))
      set(h(i,j),'fontsize',fontsize,'fontname',fontname)
      patch(X(:,xind(indij)),Y(:,yind(indij)),Z(:,indij),...
                    'EdgeColor',edgecolor,'FaceColor',facecolortype);
      hold on
        plot(h(i,j),abx([1 2 2 1 1]),aby([1 1 2 2 1]),'k-')
      hold off
      set(h(i,j),'xlim',                   abx,...
                 'ylim',                   aby,...
                 'ticklength',             [0 0], ...
                 'clim',                   clim,...
                 'ActivePositionProperty', 'OuterPosition');
      if squareplot, set(h(i,j),'DataAspectRatio',[1 1 1]); end
      if j>1
        set(h(i,j),'yticklabel',' '); 
      else
        %ylabel(labels{2})
      end
      if i>1, 
        set(h(i,j),'xticklabel',' '); 
      else
        %xlabel(labels{1})
      end
      box on
    end
  end  
  
  
  % Put  ticks, xlabels at bottom of figure and w labels at top
  for j=1:ncols
    if ~noticklabels
       set(h(1,j),'xcolor',[0 0 0]);
    end
    pos=get(h(1,j),'position');
    if nrows>1 && ncols>1
      axes('parent',hpanels(3), 'Position', [pos(1) 0 pos(3) 1],...
         'xcolor',[1 1 1], 'ycolor',[1 1 1],'xtick',[],'ytick',[],...
         'ticklength',[0 0]);
      text(.5,1,labels{1},'HorizontalAlignment','center','VerticalAlignment','top');
    else
      set(gcf,'CurrentAxes',h(1,j))
      xlabel(labels{1})
    end
    if zlabels
      axes('parent',hpanels(5), 'Position', [pos(1) 0 pos(3) 1],...
           'xcolor',[1 1 1], 'ycolor',[1 1 1],'xtick',[],'ytick',[],'xticklabel','','yticklabel','');
      if all(round(Sc)==Sc)
        if size(Sc,2)>1
          lab=[labels{3} ' = [' sprintf('%1i,',Sc(j,1:end-1)) sprintf('%1i],',Sc(j,end))]; 
        else
          lab=[labels{3} ' = ' sprintf(' %1i',Sc(j,:))]; 
        end
      else
        if size(Sc,2)>1
          lab=[labels{3} ' = [' sprintf('%1.2f,',Sc(j,1:end-1)) sprintf('%1.2f]',Sc(j,end))]; 
        else
          lab=[labels{3} ' = ' sprintf(' %1.2f',Sc(j,:))]; 
        end
      end
      text(.5,0,lab,'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
  end
  
  
  % Put  ticks, ylabels and w labels on left edge of figure
  for i=1:nrows
    if ~noticklabels
       set(h(i,1),'ycolor',[0 0 0]);
    end
    pos=get(h(i,1),'position');
    if nrows>1 && ncols>1
      axes('parent',hpanels(4), 'Position', [0 pos(2) 1 pos(4)],...
         'xcolor',[1 1 1], 'ycolor',[1 1 1],'xtick',[],'ytick',[],...
         'ticklength',[0 0]);
      text(1,.5,labels{2},'HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',90);
    else
      set(gcf,'CurrentAxes',h(i,1))
      ylabel(labels{2})
    end
    if wlabels
      axes('parent',hpanels(6), 'Position', [0 pos(2) 1 pos(4)],...
           'xcolor',[1 1 1], 'ycolor',[1 1 1],'xtick',[],'ytick',[],'ytick',[],'xticklabel','','yticklabel','')
      if all(round(Sr)==Sr)
        if size(Sr,2)>1
          lab=[labels{4} ' = [' sprintf(' %1i',Sr(i,:))  sprintf(']')];
        else
          lab=[labels{4} ' = ' sprintf(' %1i',Sr(i,:))];
        end
      else
        if size(Sr,2)>1
          lab=[labels{4} ' = [' sprintf('%1.2f,',Sr(i,1:end-1)) sprintf('%1.2f]',Sr(i,end))];
        else
          lab=[labels{4} ' = ' sprintf(' %1.2f',Sr(i,:))];
        end
      end
      text(0.5,.5,lab,'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90);
    end
  end
  

  
  if addlegend 
    if colorbartype
      if vertical
        hc=vcolorbar(hpanels(2),clim,legendtitle);
      else
        hc=hcolorbar(hpanels(1),clim,legendtitle);
      end
    else
      if vertical
        hc=vlegend(hpanels(2),vvals,clim,legendlabels,legendtitle);
      else
        hc=hlegend(hpanels(1),vvals,clim,legendlabels,legendtitle);
      end
    end
  end
  set(gcf,'CurrentAxes',h(1,1))
  
function hl=vlegend(h,vvals,clim,labels,legendtitle)
  hl=axes('parent',h);
  nv=length(labels);
  if nv==0
    labels=cell(clim(2)-clim(1)+1,1); 
    nv=length(labels);
    for i=1:nv; labels{i}=num2str(clim(1)+i-1); end
  end
  if length(vvals)==nv
    ii=1:nv;
    u=(clim(2)-clim(1))/(nv-1)/4;
    X=[-u;-u;u;u]*ones(1,nv);
    Y=[ii-u;ii+u;ii+u;ii-u];
    Z=ones(4,1)*linspace(clim(1),clim(2),nv);
    boxon=1;
    if boxon, boxcolor=[0 0 0]; set(hl,'box','on');
    else      boxcolor=[1 1 1]; set(hl,'box','off'); %#ok<UNRCH>
    end
    set(hl,'ylim',                   [1-3*u,nv+3*u],...
           'box',                    'on',...
           'xcolor',                 boxcolor,...
           'ycolor',                 boxcolor,...
           'ticklength',             [0 0], ...
           'xticklabel',             ' ',...
           'yticklabel',             ' ',...
           'clim',                   clim,...
           'DataAspectRatio',        [1 1 1],...
           'position',               [0.05 0.05 .9 .9],...
           'ActivePositionProperty', 'OuterPosition');
    patch(X,Y,Z,'edgecolor',[0 0 0]);
    maxx=0; maxi=0;
    for i=1:nv
      %text(0,i-0.5,labels{i},'VerticalAlignment','middle','HorizontalAlignment','center')
      ht=text(2*u,i,labels{i},'VerticalAlignment','middle','HorizontalAlignment','left');
      ex=get(ht,'extent');
      if maxx<ex(3), maxx=ex(3); maxxi=ht; end
    end
    % shouldn't need to do this!
    for i=1:5
      ex=get(maxxi,'extent');
      set(hl,'xlim',[-3*u,ex(1)+ex(3)+3*u])
    end
    if ~isempty(legendtitle), title(legendtitle); end
  end

function hl=hlegend(h,vvals,clim,labels,legendtitle)
  hl=axes('parent',h);
  nv=length(labels);
  if nv==0
    labels=cell(clim(2)-clim(1)+1,1); 
    nv=length(labels);
    for i=1:nv; labels{i}=num2str(clim(1)+i-1); end
  end
  if length(vvals)<=nv
    ii=1:nv;
    u=(clim(2)-clim(1))/(nv-1)/4;
    Y=[-u;-u;u;u]*ones(1,nv);
    X=[ii-u;ii+u;ii+u;ii-u];
    Z=ones(4,1)*linspace(clim(1),clim(2),nv);
    boxon=1;
    if boxon, boxcolor=[0 0 0]; set(hl,'box','on');
    else      boxcolor=[1 1 1]; set(hl,'box','off'); %#ok<UNRCH>
    end
    set(hl,'xlim',                   [1-3*u,nv+3*u],...
           'ylim',                   [-4*u,2*u],...
           'box',                    'on',...
           'xcolor',                 boxcolor,...
           'ycolor',                 boxcolor,...
           'ticklength',             [0 0], ...
           'xticklabel',             ' ',...
           'yticklabel',             ' ',...
           'clim',                   clim,...
           'DataAspectRatio',        [1 1 1],...
           'Position',               [0.05 0.1 .95 .9],...
           'ActivePositionProperty', 'OuterPosition');
    patch(X,Y,Z,'edgecolor',[0 0 0]);
    maxx=0;
    for i=1:nv
      ht=text(i,-2.5*u,labels{i},'VerticalAlignment','middle','HorizontalAlignment','center');
      ex=get(ht,'extent');
      maxx=max(maxx,ex(3));
    end
    if ~isempty(legendtitle), title(legendtitle); end
  end

function hl=vcolorbar(h,clim,legendtitle)
  hl=axes('parent',h);
  C=colormap;
  nv=size(C,1);
  ii=linspace(clim(1),clim(2),nv);
    u=(clim(2)-clim(1))/(nv-1)/2;
    X=[-u;-u;u;u]*ones(1,nv);
    Y=[ii-u;ii+u;ii+u;ii-u];
    Z=ones(4,1)*ii;
    set(hl,'ylim',                   [clim(1)-u,clim(2)+u],...
           'box',                    'on',...
           'xcolor',                 [1 1 1],...
           'ycolor',                 [0 0 0],...
           'YAxisLocation',          'right',...
           'ticklength',             [0 0], ...
           'clim',                   clim,...
           'ActivePositionProperty', 'OuterPosition');
    patch(X,Y,Z,'edgecolor','none');
    pos=get(hl,'position'); pos(1)=0.3; pos(3)=0.15;
    set(hl,'xlim',[-u,u],'position',pos)
    title(legendtitle)

    
  function hl=hcolorbar(h,clim,legendtitle)
  hl=axes('parent',h);
  C=colormap;
  nv=size(C,1);
  ii=linspace(clim(1),clim(2),nv);
    u=(clim(2)-clim(1))/(nv-1)/2;
    Y=[-u;-u;u;u]*ones(1,nv);
    X=[ii-u;ii+u;ii+u;ii-u];
    Z=ones(4,1)*ii;
    set(hl,'xlim',                   [clim(1)-u,clim(2)+u],...
           'box',                    'on',...
           'ycolor',                 [1 1 1],...
           'xcolor',                 [0 0 0],...
           'ticklength',             [0 0], ...
           'clim',                   clim,...
           'ActivePositionProperty', 'OuterPosition');
    patch(X,Y,Z,'edgecolor','none');
    if length(legendtitle)>0
      h=title(legendtitle);
      pos=get(hl,'position'); pos(2)=0.6; pos(4)=0.15;
    else
      pos=get(hl,'position'); pos(2)=0.75; pos(4)=0.2;
    end
    set(hl,'ylim',[-u,u],'position',pos)