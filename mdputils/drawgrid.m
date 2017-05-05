% drawgrid Draws a grid on the current axis
% USAGE
%   drawgrid(x,y,sym);
% INPUTS
%   x   : vector of x values
%   y   : vector of y values
%   sym : (optional) symbol used to control color and line type [default: 'k:']

function drawgrid(x,y,sym)
if nargin<3, sym='k:'; end
%xx=get(gca,'xlim');
%yy=get(gca,'ylim');
xx=[min(x) max(x)];
yy=[min(y) max(y)];
hold on
for i=1:length(x)
  plot([x(i) x(i)],yy,sym)
end
for i=1:length(y)
  plot(xx,[y(i) y(i)],sym)
end

xx=get(gca,'xlim');
yy=get(gca,'ylim');
for i=1:length(x)
  plot([x(i) x(i)],yy,sym)
end
for i=1:length(y)
  plot(xx,[y(i) y(i)],sym)
end
hold off