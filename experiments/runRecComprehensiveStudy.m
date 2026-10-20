function directory=runRecComprehensiveStudy(profile,blocks)
%RUNRECCOMPREHENSIVESTUDY 冻结协议、保存源码快照，并按独立条件块并行续跑。
%   每个文件包含全部配对试次，不将一次试次当作独立几何样本。
arguments
    profile (1,1) string = "smoke"
    blocks (1,:) string = "all"
end
root=fileparts(fileparts(mfilename('fullpath'))); cfg=config(profile);
% 事后数值复核单列批次，不用更密网格覆盖已有主实验或按表现换几何。
dense=isequal(blocks,"numerical_refinement");
if dense
    cfg.experiments.recComprehensive.numericalUpaGrids=cfg.experiments.recComprehensive.convergenceUpaGrids;
    cfg.experiments.recComprehensive.workers=cfg.experiments.recComprehensive.numericalWorkers;
end
[geometries,cases]=uav.experiment.recComprehensivePlan(cfg);
if dense
    cases=cases(string({cases.block})=="numerical" & string({cases.arrayType})=="UPA");
    for i=1:numel(cases),cases(i).id=i;cases(i).block="numerical_refinement";end
elseif ~isequal(blocks,"all")
    cases=cases(ismember(string({cases.block}),blocks));
    assert(~isempty(cases),'选择的实验区块不存在。');
    for i=1:numel(cases),cases(i).id=i;end
end
files=[dir(fullfile(root,'src','**','*.m'));dir(fullfile(root,'experiments','*.m'));dir(fullfile(root,'config.m'))];
sources=struct('path',{},'text',{}); digest=java.security.MessageDigest.getInstance('SHA-256');
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    sources(i).path=erase(string(path),string(root)+filesep);
    sources(i).text=fileread(path);
    digest.update(unicode2native(char(sources(i).path+newline+string(sources(i).text)),'UTF-8'));
end
% 派生流形含复数且可由配置重建；签名只编码原始参数，完整派生量仍保存。
parameters=rmfield(cfg,'derived');
digest.update(unicode2native(jsonencode(parameters),'UTF-8'));
digest.update(unicode2native(jsonencode(blocks),'UTF-8'));
signature=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
parent=fullfile(root,cfg.experiments.recComprehensive.outputDirectory,profile);
directory=fullfile(parent,signature(1:16)); if ~isfolder(directory),mkdir(directory);end
manifest=struct('cfg',cfg,'geometries',geometries,'cases',cases,'sources',sources, ...
    'sourceSignature',signature,'createdAt',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'matlabVersion',version,'toolboxes',ver);
manifest.selectedBlocks=blocks;
[~,commit]=system(sprintf('git -C "%s" rev-parse HEAD',root)); manifest.baseCommit=strtrim(commit);
[~,dirty]=system(sprintf('git -C "%s" status --short',root)); manifest.workingTreeStatus=dirty;
if ~isfile(fullfile(directory,'manifest.mat')),save(fullfile(directory,'manifest.mat'),'manifest','-v7.3');end
pointer="latest_run.mat";
if ~isequal(blocks,"all"),pointer="latest_"+strjoin(blocks,"_")+"_run.mat";end
save(fullfile(parent,pointer),'directory');
pending=[];
for i=1:numel(cases)
    path=fullfile(directory,sprintf('case_%06d.mat',i));
    if ~isfile(path),pending(end+1)=i;end %#ok<AGROW>
end
fprintf('协议 %s：总计%d块，待运行%d块，输出%s\n',cfg.experiments.recComprehensive.protocolVersion,numel(cases),numel(pending),directory);
if isempty(pending),return;end
pool=gcp('nocreate');
if isempty(pool),pool=parpool('Processes',cfg.experiments.recComprehensive.workers);end
pctRunOnAll maxNumCompThreads(1);
queue=parallel.pool.DataQueue; completed=0; timer=tic;
afterEach(queue,@progress);
% 重任务的耗时与轻任务差异很大，每次分配一个条件，避免尾部只剩少数进程工作。
options=parforOptions(pool,'RangePartitionMethod','fixed','SubrangeSize',1);
parfor (index=1:numel(pending),options)
    id=pending(index); condition=cases(id); geo=geometries(condition.geometryIndex);
    try
        result=uav.experiment.runRecComprehensiveCase(cfg,geo,condition);
        result.sourceSignature=signature;
        saveCase(fullfile(directory,sprintf('case_%06d.mat',id)),result);
        send(queue,struct('id',id,'ok',true));
    catch exception
        failure=struct('caseId',id,'identifier',exception.identifier,'message',exception.message, ...
            'stack',exception.stack,'sourceSignature',signature);
        saveFailure(fullfile(directory,sprintf('failure_%06d.mat',id)),failure);
        send(queue,struct('id',id,'ok',false));
    end
end
fprintf('本次运行结束，用时%.1f秒。请汇总检查缺失块、失败项及不适用项。\n',toc(timer));
    function progress(message)
        completed=completed+1;
        if ~message.ok || completed==1 || mod(completed,10)==0 || completed==numel(pending)
            fprintf('[%s] %d/%d，条件%d，成功%d，用时%.1fs\n', ...
                char(datetime('now','Format','HH:mm:ss')),completed,numel(pending),message.id,message.ok,toc(timer));
        end
    end
end

function saveCase(path,result)
%SAVECASE 临时文件写完后重命名；中途退出不会留下假完成文件。
temporary=path+".partial"; save(temporary,'result','-v7.3'); movefile(temporary,path);
end
function saveFailure(path,failure)
save(path,'failure');
end
