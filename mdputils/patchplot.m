% patchplot Generates a basic patch plot in the current figure 
% USAGE
%   h=patchplot(x,y,z,zvals,shownum);
% INPUTS
%   x       : n-vector of values from a set of evenly spaced values
%   y       : n-vector of values from a set of evenly spaced values
%   z       : n-vector of values
%   zvals   : a vector of possible values for z - 
%               this will place a number of dummy patches on the plot 
%               that are needed if a legend is added to the figure
%   shownum : places the value of z as text on the plot
% OUTPUT
%   h : handle to the patch object
% 
% Generates a set of patches with the x-data on the x-axis, y-data
% on the y-axis and uses the z-data to determine the color of each patch
% based on the current figure color map.
%
% Additional information can be added to this plot including titles, axis labels, 
% and legends. It can be also called to plot in subplot.

function h=patchplot(x,y,z,zvals,shownum)
if nargin<5 || isempty(shownum)
  shownum=false;
end
if isempty(x)
  warning('no values passed - nothing to do')
  return
end
if nargin<4, zvals=[]; end
minx=min(x); maxx=max(x);
miny=min(y); maxy=max(y);
xinc=(maxx-minx)/(numel(unique(x))-1);
yinc=(maxy-miny)/(numel(unique(y))-1);
if isnan(xinc), xinc=1; end
if isnan(yinc), yinc=1; end
n=length(x);
xx=[-xinc;xinc;xinc;-xinc]/2;
yy=[-yinc;-yinc;yinc;yinc]/2;

if n==3
  x=[x(1);x];
  y=[y(1);y];
  z=[z(1);z];
  n=4;
end

if isempty(zvals)
  minz=min(z);
  maxz=max(z);
  clim=[minz maxz];
else
  for i=1:length(zvals)
    x1=minx-[i i+1]-1; y1=miny-[i i+1]-1;
    patch(xx*x1+ones(4,1)*x1,yy*y1+ones(4,1)*y1,[zvals(i) zvals(i)], ...
        'Edgecolor','none','linestyle','none');
  end
  clim=[min(zvals) max(zvals)];
end

ch=patch(xx*ones(1,n)+ones(4,1)*x',yy*ones(1,n)+ones(4,1)*y',z', ...
    'Edgecolor','none','linestyle','none');
if ~isnan(xinc)
  xlim([minx-xinc/2 maxx+xinc/2])
end
if ~isnan(yinc)
  ylim([miny-yinc/2 maxy+yinc/2])
end
if diff(clim)>0, set(gca,'clim',clim); end

if shownum
  for i=1:length(z)
    text(x(i),y(i),num2str(z(i)),'HorizontalAlignment','center');
  end
end

if nargout>0, h=ch; end
