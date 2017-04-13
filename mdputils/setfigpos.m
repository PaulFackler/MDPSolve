% prints out a code string top set the current figure position to it's current place
set(gcf,'units','normalized')
pos=get(gcf,'position');
fprintf('set(gcf,''units'',''normalized'',''position'',[%5.3f %5.3f %5.3f %5.3f])\n',pos)