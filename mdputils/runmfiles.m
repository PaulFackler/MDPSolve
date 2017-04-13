% runmfiles Runs all of the m-files in the current directory
set(0,'DefaultFigureUnits','normalized','DefaultFigurePosition',[.5 .25 .5 .5])
files = dir;
k=0;
names=cell(size(files,1),1);
for i=1:size(files)
  fn=files(i).name; 
  if length(fn)>2 && strcmp(fn(end-1:end),'.m')
    if ~strcmp(fn(1:end-2),mfilename)
      k=k+1';
      names{k}=fn(1:end-2);
    end
  end
end
names(k+1:end)=[];

for dispi=1:k
  clc
  close all
  disp(['Running: ' names{dispi}])
  clearvars -except names dispi k 
  try
     eval(names{dispi})
  catch %#ok<CTCH>
     disp('error encountered - continuing to next file')
  end
  disp(' ')
  disp(['Finished running: ' names{dispi} ' -- Press any key to continue'])
  ii=get(0,'Children');
  if ~isempty(ii), figure(min(ii)); end
  pause
end
   
clc