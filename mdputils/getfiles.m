% if profile is on this lists all non-Matlab files run since
% profiler was turned on.
function getfiles
stats = profile('info');
x=stats.FunctionTable;
y=cell(1,1);
k=1;
for i=1:size(x,1)
  s=x(i).CompleteName;
  if isempty(strfind(s,'C:\Program Files\MATLAB\')) && ...
     isempty(strfind(s,'com.mathworks')) && ...
     isempty(strfind(s,'java.io.File')) && ...
     isempty(strfind(s,'@('))
    y{k}=s;
    k=k+1;
  end
end
y=sort(y);
for i=1:numel(y)
  disp(y{i})
end