function startup()
%STARTUP Add only repository runtime paths required by this simulator.

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'src'));
addpath(fullfile(rootDir, 'experiments'));
addpath(fullfile(rootDir, 'tests'));
end
