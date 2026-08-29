function result = run_all(profile, options)
%RUN_ALL 按固定顺序运行全部测试和五组实验。
%   result = run_all(profile)
%   result = run_all(profile, RunTests=false, RunExperiments=true)
%
%   输出result保存profile、测试结果、实验名称和每组实验的产物路径。
%   默认先测试后实验：任何测试失败都会由assertSuccess立即终止，避免在
%   已知代码不满足契约时继续生成看似正式的结果。

arguments
    profile (1,1) string = "preliminary"
    options.RunTests (1,1) logical = true
    options.RunExperiments (1,1) logical = true
end

startup;
cfg = config(profile);
% run_all属于批处理入口，每组实验完成导出后立即关闭图窗。
cfg.plot.closeAfterExport = true;
rootDir = cfg.output.rootDir;
result = struct('profile', cfg.profile, 'tests', [], ...
    'experiments', strings(0, 1), 'artifacts', {cell(0, 1)});

if options.RunTests
    % IncludeSubfolders保证五个阶段的测试类全部被发现。
    result.tests = runtests(fullfile(rootDir, 'tests'), 'IncludeSubfolders', true);
    assertSuccess(result.tests);
end

if options.RunExperiments
    % entryPoints提供唯一的实验顺序，防止文件系统排序改变执行顺序。
    entryPoints = uav.experiment.entryPoints();
    for index = 1:numel(entryPoints)
        experimentResult = feval(entryPoints(index), cfg);
        result.experiments(end + 1, 1) = entryPoints(index);
        result.artifacts{end + 1, 1} = experimentResult.artifacts;
    end
end
end
