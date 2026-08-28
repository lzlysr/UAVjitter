function ensureDirectories(output)
%ENSUREDIRECTORIES 创建配置指定的缓存与结果目录。
%   已存在的目录保持不变；本函数不删除、不清空任何用户文件。

fields = {'rawDir', 'tableDir', 'figureDir', 'logDir', 'cacheDir'};
for index = 1:numel(fields)
    pathValue = output.(fields{index});
    if ~isfolder(pathValue)
        mkdir(pathValue);
    end
end
end
