function resumeRecComprehensiveStudy(directory)
%RESUMERECCOMPREHENSIVESTUDY 从冻结清单续跑，一条件一任务，避免重任务尾部空等。
%   不重新读取当前config；只使用原清单和未变化的计算核心。
load(fullfile(directory,'manifest.mat'),'manifest');
root=fileparts(fileparts(mfilename('fullpath')));
for i=1:numel(manifest.sources)
    s=manifest.sources(i); normalized=replace(string(s.path),"\","/");
    if startsWith(normalized,"src/") && ~endsWith(normalized,"recComprehensivePlan.m")
        assert(strcmp(fileread(fullfile(root,s.path)),s.text), ...
            'uav:experiment:ChangedCore','计算核心与冻结版本不同，拒绝混用结果。');
    end
end
pending=find(~arrayfun(@(id)isfile(fullfile(directory,sprintf('case_%06d.mat',id))),1:numel(manifest.cases)));
resumeMetadata=struct('pendingCases',pending,'startedAt',char(datetime('now')), ...
    'sourceSignature',manifest.sourceSignature,'resumeSource',fileread(mfilename('fullpath')+".m"));
save(fullfile(directory,'resume_metadata.mat'),'resumeMetadata');
if isempty(pending),return;end
pool=gcp('nocreate');if isempty(pool),pool=parpool('Processes',manifest.cfg.experiments.recComprehensive.workers);end
pctRunOnAll maxNumCompThreads(1);
futures=parallel.FevalFuture.empty;
for i=1:numel(pending)
    c=manifest.cases(pending(i));geo=manifest.geometries(c.geometryIndex);
    futures(i)=parfeval(pool,@uav.experiment.runRecComprehensiveCase,1,manifest.cfg,geo,c);
end
timer=tic;
for completed=1:numel(pending)
    [index,result]=fetchNext(futures);id=pending(index);
    result.sourceSignature=manifest.sourceSignature;
    path=fullfile(directory,sprintf('case_%06d.mat',id));temporary=path+".partial";
    assert(~isfile(path),'条件已由其他进程完成；禁止覆盖并行写入。');
    save(temporary,'result','-v7.3');movefile(temporary,path);
    if completed==1||mod(completed,5)==0||completed==numel(pending)
        fprintf('续跑%d/%d，完成条件%d，用时%.1fs\n',completed,numel(pending),id,toc(timer));
    end
end
end
