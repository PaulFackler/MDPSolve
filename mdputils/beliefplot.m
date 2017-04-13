% beliefplot Generates a basic patch plot in the current figure 
% USAGE
%   h=beliefplot(b,z,zvals,tsymbol,loffset,toffset);
% INPUTS
%   b       : nx3 matrix of belief weights (from simplexgrid)
%   z       : n-vector of values
%   zvals   : a vector of possible values for z - this will place a number of dummy
%             patches on the plot that are needed if a legend is added to the figure
%   tsymbol : symbol for axis labeling (e.g., 'b' or 'w')
%   loffset : distance from axis to axis labels [default: 0.125]
%   toffset : distance from axis to tick labels [default: 0.09]
% OUTPUT
%   h : handle to the patch object
% 
% Generates a set of patches with the x-data on the x-axis, y-data
% on the y-axis and uses the z-data to determine the color of each patch
% based on the current figure color map.
%
% Additional information can be added to this plot including titles, axis labels, 
% and legends. It can be also called to plot in subplot.

function h=beliefplot(b,z,zvals,tsymbol,loffset,toffset)
if nargin<6 || isempty(toffset)
  toffset=0.09;
end
if nargin<5 || isempty(loffset)
  loffset=0.125;
end
if nargin<4 
  tsymbol='b';
end
if nargin<3, zvals=[]; end

fontname='Times New Roman';

cc=sqrt(3)/4;   % used for equilateral triangle calculations
x=b*[0;1;0.5];  y=b*[0;0;cc];
minx=min(x); maxx=max(x);
miny=min(y); maxy=max(y);
xinc=4*(maxx-minx)/(numel(unique(x))-1);
yinc=(maxy-miny)/(numel(unique(y))-1);
n=length(x);
xx=[-xinc;xinc;xinc;-xinc]/2;
yy=[-yinc;-yinc;yinc;yinc]/2;

% determine the range of the plotted values
if isempty(zvals)
  minz=min(z);
  maxz=max(z);
  clim=[minz maxz];
else
  % add dummy patches to enable legends
  for i=1:length(zvals)
    x1=minx-[i i+1]-1; y1=miny-[i i+1]-1;
    h=patch(xx*x1+ones(4,1)*x1,yy*y1+ones(4,1)*y1,[zvals(i) zvals(i)], ...
        'Edgecolor','none','linestyle','none');
  end
  clim=[min(zvals) max(zvals)];
end

% plot squares around each belief point
h=patch(xx*ones(1,n)+ones(4,1)*x',yy*ones(1,n)+ones(4,1)*y',z', ...
    'Edgecolor','none','linestyle','none');
% defines the relative color 
if diff(clim)>0, set(gca,'clim',clim); end
% is this needed ?
xlim([minx-xinc/2 maxx+xinc/2])
ylim([miny-yinc/2 maxy+yinc/2])

bc=get(gcf,'color');
set(gca,'xcolor',bc,'ycolor',bc)
set(gca,'DataAspectRatio', [2 1 1])
set(gca,'xtick',[],'ytick',[])

sloffset=loffset/sqrt(2);
stoffset=toffset/sqrt(2);
text(0.5,-cc*loffset/2,[tsymbol '_2'],'HorizontalAlignment','center','FontName',fontname);
text(0,-cc*toffset/2,'0','HorizontalAlignment','center','FontName',fontname);
text(1,-cc*toffset/2,'1','HorizontalAlignment','center','FontName',fontname);

text(0.25-sloffset,cc/2+cc*sloffset/2,[tsymbol '_1'],'rotation',60,'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname);
text(0.5-stoffset,cc+cc*stoffset/2,'0','rotation',60,'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname);
text(-stoffset,cc*stoffset/2,'1','rotation',60,'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname);

text(0.75+sloffset,cc/2+cc*sloffset/2,[tsymbol '_3'],'rotation',300,'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname);
text(0.5+stoffset,cc+cc*stoffset/2,'1','rotation',300,'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname);
text(1+stoffset,cc*stoffset/2,'0','rotation',300,'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname);

return
hold on
%plot([0 1 0.5 0],[0 0 cc 0],'k')
%plot([.25 0.5 .75 .25],[.5 0 .5 .5]*cc,'k')
hold off

%axis([0 1 0 cc])
%axis square
