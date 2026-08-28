function startup()
%STARTUP 将仿真系统所需目录加入MATLAB搜索路径。
%   这里只做路径管理，不创建配置、不运行试验，也不改变随机数状态。
%   +uav 是MATLAB包目录，因此只需加入其父目录src。

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'src'));
addpath(fullfile(rootDir, 'experiments'));
addpath(fullfile(rootDir, 'tests'));
end
