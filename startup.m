% Adds the MDPSOLVE path so MATLAB can find the toolbox and demos.
% Other things can be added to this file to set the defaults you like.
% Putting this file in the right place allows it to be run automatically
%   when you start MATLAB. See MATLAB documentation on startup file.

% change the path to reflect where your version of MDPSolve is located
addpath(genpath('C:\MDPSOLVE'))

format compact 

% default figure settings
set(0,'DefaultFigureColor',[1 1 1])
set(0,'DefaultSurfaceFaceColor','interp','DefaultSurfaceEdgeColor','interp')
set(0,'defaultAxesFontName', 'times new roman')
set(0,'defaultTextFontName', 'times new roman')
set(0,'defaultAxesFontSize', 12)
set(0,'defaultTextFontSize', 12)
set(0, 'DefaultFigureUnits', 'inches');
set(0, 'DefaultFigurePaperPositionMode', 'auto');
set(0, 'DefaultFigurePosition', [7 0 8 6]);
