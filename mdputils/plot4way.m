function plot4way(X,Y,xvals,yvals,xlabs,titlestring,legend)

if iscell(X)
  X1=X{1};
  X2=X{2};
  X3=X{3};
  X4=X{4};
else
  X1=X(:,1);
  X2=X(:,2);
  X3=X(:,3);
  X4=X(:,4);
end
if nargin<3 || isempty(xvals)
  xvals{1}=[]; 
  xvals{2}=[]; 
  xvals{3}=[]; 
  xvals{4}=[]; 
end
if nargin<6,  titlestring=''; end
if isempty(xvals{3}), xvals{3}=unique(X3); end
if isempty(xvals{4}), xvals{4}=unique(X4); end
if nargin<4, yvals=[]; end
if nargin<5 || isempty(xlabs) 
  xlabs{1}='X1'; 
  xlabs{2}='X2'; 
  xlabs{3}='X3=%1.1f'; 
  xlabs{4}='X4=%1.1f'; 
end

xstrip=0.0125;
ystrip=0.0125;
n=length(xvals{3});
m=length(xvals{4});
xlen=0.8/n;
ylen=0.8/m;
h=cell(m,n);
set(gcf,'units','normalized')
for i=1:n
  ind3=X3==xvals{3}(i);
  for j=1:m
    ind=ind3 & X4==xvals{4}(j);
    h{i,j}=axes('position',[0.1+(i-1)*xlen 0.1+(j-1)*ylen xlen-xstrip ylen-ystrip]);
    patchplot(X1(ind),X2(ind),Y(ind),yvals);
    if j>1, set(gca,'XTick',[]); 
    else
      xlabel(xlabs{1});
      if ~isempty(xvals{1})
        set(gca,'Xtick',xvals{1});
      end
    end
    if i>1, set(gca,'YTick',[]); 
    else
      ylabel(xlabs{2});
      if ~isempty(xvals{2})
        set(gca,'Ytick',xvals{2});
      end
      axes('position',[0 0.1+(j-1)*ylen 0.1 ylen],'visible','off');
      hh=text(0.25,0.5,sprintf(xlabs{4},xvals{4}(j)),...
        'VerticalAlignment','middle','Rotation',90,'HorizontalAlignment','center');
    end
  end
  axes('position',[0.1+(i-1)*xlen 0 xlen 0.1],'visible','off');
  hh=text(0.5,0.25,sprintf(xlabs{3},xvals{3}(i)),...
        'VerticalAlignment','middle','Rotation',0,'HorizontalAlignment','center');
end
axes(h{n,m})
if legend<0
  hh=colorbar;
  pos=get(hh,'outerposition');
  set(hh,'outerposition',[1-2*pos(3) 0.25 pos(3) 0.5])
end
if ~isempty(titlestring)
  axes('position',[0.1 .9 .8 0.1],'visible','off');
  hh=text(0.5,0.25,titlestring,...
        'VerticalAlignment','middle','Rotation',0,'HorizontalAlignment','center');
  fs=get(hh,'FontSize');
  set(hh,'FontSize',2*fs)
end