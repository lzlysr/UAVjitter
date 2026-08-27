function result = run_all(profile, options)
%RUN_ALL Run the test suite and each available experiment entry point in order.

arguments
    profile (1,1) string = "preliminary"
    options.RunTests (1,1) logical = true
    options.RunExperiments (1,1) logical = true
end

startup;
cfg = config(profile);
rootDir = cfg.output.rootDir;
result = struct('profile', cfg.profile, 'tests', [], 'experiments', strings(0, 1));

if options.RunTests
    result.tests = runtests(fullfile(rootDir, 'tests'), 'IncludeSubfolders', true);
    assertSuccess(result.tests);
end

if options.RunExperiments
    experimentFiles = dir(fullfile(rootDir, 'experiments', 'experiment_*.m'));
    [~, order] = sort({experimentFiles.name});
    experimentFiles = experimentFiles(order);
    for index = 1:numel(experimentFiles)
        [~, entryPoint] = fileparts(experimentFiles(index).name);
        feval(entryPoint, cfg);
        result.experiments(end + 1, 1) = string(entryPoint); %#ok<AGROW>
    end
end
end
