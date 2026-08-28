function ensureDirectories(output)
%ENSUREDIRECTORIES Create configured generated-output directories.

fields = {'rawDir', 'tableDir', 'figureDir', 'logDir', 'cacheDir'};
for index = 1:numel(fields)
    pathValue = output.(fields{index});
    if ~isfolder(pathValue)
        mkdir(pathValue);
    end
end
end
