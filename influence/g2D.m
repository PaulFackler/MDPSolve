function D=g2D(g,s,x,e,xelist,options)
print=0;
order=[];
if nargin<6, options=[]; end
if isnumeric(options)
  print=options;
else
  getopts(options, ...
         'print',       0, ...
         'order',       []); 
end
g=inputcheck(g,'g','function_handle');
s=inputcheck(s,'s','numeric');
x=inputcheck(x,'x','numeric');
e=inputcheck(e,'e','struct');
xelist=inputcheck(xelist,'xlist','numeric');

ds=length(s);
dx=length(x);
de=length(e);

xlist=cell(1,ds);
elist=cell(1,ds);
for i=1:ds
  if length(xelist{i})~=nargin(g{i})
    error(['xelist has a different number of inputs than g for variable ' num2str(i)])
  end
  elist{i}=-xelist{i}(xelist{i}<0);
  if any(elist{i}>de)
    error(['xelist{' num2str(i) '} contains a value greater than the # of noise variables'])
  end
  xlist{i}=xelist{i}(xelist{i}>0);
  if any(xlist{i}>dx)
    error(['xelist{' num2str(i) '} contains a value greater than the # of X variables'])
  end
end
names=cell(1,dx+de);
cpdnames=cell(1,dx+de);
for i=1:ds,    names{i}   =['S' num2str(i)]; cpdnames{i}      =['x{' num2str(i) '}'];    end
for i=1:dx-ds, names{ds+i}=['A' num2str(i)]; cpdnames{ds+i}   =['x{' num2str(ds+i) '}']; end
for i=1:de,    names{dx+i}=['e' num2str(i)]; cpdnames{dx+i}   =['e{' num2str(i) '}'];    end
for i=1:ds,    cpdnames{dx+de+i}   =['g{' num2str(i) '}'];    end
D=[];
for i=1:ds
  D=add2diagram(D,names{i},'s',1,{},s{1});
end
for i=1:dx-ds
  D=add2diagram(D,names{ds+i},'a',1,{},x{ds+i});
end
for i=1:de
  D=add2diagram(D,names{dx+i},'c',1,{},e{i});
end
for i=1:ds
  di=length(xelist{i});
  parentsi=cell(1,di);
  for j=1:di
    if xelist{i}(j)>0, parentsi{j}=names{xelist{i}(j)};
    else               parentsi{j}=names{dx-xelist{i}(j)};
    end
  end
  D=add2diagram(D,[names{i} '+'],'f',1,parentsi,g{i});
end
if length(g)>ds
  i=ds+1;
  di=length(xelist{i});
  parentsi=cell(1,di);
  for j=1:di
    if xelist{i}(j)>0, parentsi{j}=names{xelist{i}(j)};
    else               parentsi{j}=names{dx-xelist{i}(j)};
    end
  end
  D=add2diagram(D,['U'],'u',1,parentsi,g{i});
  cpdnames{dx+de+ds+1}   =['g{' num2str(i) '}'];
end
D.cpdnames=cpdnames;
if print
  printD(D,1);
end

function x=inputcheck(x,xname,class)
if ~iscell(x)
  if isa(x,'class')
    x={x};
  else
    error([xname ' must be a single element or a cell array of elements of class ' class])
  end
end