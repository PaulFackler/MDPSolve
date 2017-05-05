% patchplot Generates a basic patch plot in the current figure 
% USAGE
%   h=patchplot(x,y,z,zvals,shownum);
% INPUTS
%   x : n-vector of values from a set of evenly spaced values
%   y : n-vector of values from a set of evenly spaced values
%   z : n-vector of values
%   zvals : a vector of possible values for z - this will place a number of dummy
%           patches on the plot that are needed if a legend is added to the figure
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

[xi,junk,indx]=unique(x);
[yi,junk,indy]=unique(y);

minx=min(xi); maxx=max(xi);
miny=min(yi); maxy=max(yi);

indx=indx+1;
indy=indy+1;
xi=[2*xi(1)-xi(2); xi; 2*xi(end)-xi(end-1)];
yi=[2*yi(1)-yi(2); yi; 2*yi(end)-yi(end-1)];


if nargin<4, zvals=[]; end
if isempty(zvals)
  minz=min(z);
  maxz=max(z);
  clim=[minz maxz];
else
  for i=1:length(zvals)
    patch([minx;maxx;maxx;minx],[miny;miny;maxy;maxy],ones(4,1)*zvals(i), ...
        'Edgecolor','none','linestyle','none');
  end
  clim=[min(zvals) max(zvals)];
end

xx=[(xi(indx)+xi(indx-1))'/2;
    (xi(indx)+xi(indx+1))'/2;
    (xi(indx)+xi(indx+1))'/2;
    (xi(indx)+xi(indx-1))'/2];
yy=[(yi(indy)+yi(indy-1))'/2;
    (yi(indy)+yi(indy-1))'/2;
    (yi(indy)+yi(indy+1))'/2;
    (yi(indy)+yi(indy+1))'/2];
h=patch(xx,yy,z', ...
    'Edgecolor','none','linestyle','none');
xlim([(xi(1)+xi(2))/2 (xi(end)+xi(end-1))/2])
ylim([(yi(1)+yi(2))/2 (yi(end)+yi(end-1))/2])
if diff(clim)>0, set(gca,'clim',clim); end

if shownum
  for i=1:length(z)
    text(x(i),y(i),num2str(z(i)),'HorizontalAlignment','center');
  end
end

if nargout==0, h=[]; end
