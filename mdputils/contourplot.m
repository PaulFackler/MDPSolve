% cotourplot Generates a basic contour plot in the current figure 
% USAGE
%   h=contourplot(x,y,z,zvals,shownum);
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

function [varargout]=contourplot(x,y,z,zvals,shownum)
if nargin<5 || isempty(shownum)
  shownum=false;
end
if isempty(x)
  warning('no values passed - nothing to do')
  return
end
if nargin<4, zvals=[]; end
if isempty(zvals)
  minz=min(z);
  maxz=max(z);
  zvals=(maxz+minz)/2; % display median value
end

for i=1:length(zvals)
  zi=zvals(i);
  ii=find(z==zi);
  plot(x(ii),y(ii),'k-')
  hold on
end

if shownum
end
hold off
if nargout>0, varargout{1}=gca; end

% workaround for bug in OpenGL renderer
%ct = camtarget;
%camtarget([ct(1)+0.001*ct(1) ct(2)+0.001 ct(3)]);


set(gca,'ZTick',[]);
