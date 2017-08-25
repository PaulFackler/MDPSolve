% MDPMEXALL Creates MEX files for MDPSOLVE toolbox
% Run this function once after the toolbox is installed
% Requires a C compiler
% Type help mex for more information on creating mex files
% You may want to set compiler options to change speed/memory 
%   default settings
%
% This program need be run only once
%
% Note: this will compile all of the C files in the 
%   MDPSOLVE disrectory and its subdirectories

function mdpmexall

% save the name of the current default directory
currentdir=cd;
% switch to the MDPSOLVE directory
mdpdir=which('mdpmexall');
mdpdir=stripdir(mdpdir);
% process subdirectories
cd(mdpdir)
processdir
% switch back to original default directory
cd(currentdir)

function mdpdir=stripdir(mdpdir)
while mdpdir(end)~='\' && mdpdir(end)~='/', mdpdir(end)=[]; end


function processdir
currentdir=cd;
fn=dir;
for i=3:length(fn)
  if fn(i).isdir
    if ~strcmp(fn(i).name(1),'@') && ...
      ~strcmp(fn(i).name,'kdtree') && ...
      ~strcmp(fn(i).name,'tprod')
      if(isunix)
          cd(['./' fn(i).name])
      else
           cd(['.\' fn(i).name])
      end
      processdir
      cd(currentdir)
    end
  elseif strcmp(fn(i).name(end-1:end),'.c')
    % mex all C files in the mdputils subdirectory
    for ii=1:10
      try
        eval(['mex -largeArrayDims -lmwblas ' fn(i).name])
        break
      catch errobj
        if strfind(errobj.message,...
         'mt : general error c101008d: Failed to write the updated manifest to the resource of file')>0
        else
          rethrow(errobj);
        end
      end
    end
    if ii>=10, disp(['Manifest file not written of ' fn(i).name]); end
    disp(['mex file created for ' cd '\' fn(i).name])
  end
end
