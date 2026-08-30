# UPA二维INCM诊断与协方差拟合修正实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不修改原始Gu–Leshem Rec、不读取真值的前提下，诊断其二维UPA直接迁移失败机理，并实现同时支持8和64快拍的 `upa_covfit_incm`。

**Architecture:** 实验9分成“零抖动网格/测度诊断”和“0°/4°/8°抖动验证”两条配对数据流。原始Rec只在样本协方差可逆时运行；新方法以加载Capon谱提取期望区与干扰区候选方向，通过非负Hermitian协方差拟合估计候选功率，并把干扰候选与独立噪声项重构为INCM。

**Tech Stack:** MATLAB R2024a、`matlab.unittest`、Optimization Toolbox的 `lsqnonneg`、现有 `+uav` 包、确定性Threefry子流与结果事务导出器。

**Spec:** `docs/superpowers/specs/2026-08-30-upa-incm-diagnosis-design.md`

## Global Constraints

- 所有新文档和代码注释使用中文。
- 原始 `src/+uav/+incm/guLeshem2012Rec.m` 不得修改。
- `L=8<M=9` 时原始Rec必须记录 `singular-not-applicable`，不得调用伪逆或加载后冒用原名称。
- `upa_covfit_incm`不得读取 `obs.truth`、真实干扰数、真实方向或真实姿态。
- 所有新增数值参数只能进入 `cfg.literature.guLeshem2012.upaDiagnosis`。
- 同一试次的8快拍观测必须逐元素等于64快拍观测的前8列。
- `smoke`只验证链路；性能判断使用 `preliminary`配对蒙特卡洛。
- 不能通过削弱ESP/WCO、改变随机种子或只报告有利点满足性能目标。

---

### Task 1: 双快拍配置合同与三类UPA求积网格

**Files:**
- Modify: `config.m`
- Modify: `src/+uav/+config/validateGuLeshem2012.m`
- Create: `src/+uav/+spectrum/upaIntegrationGrid.m`
- Create: `tests/TestTask8UpaIncmDiagnosis.m`

**Interfaces:**
- Consumes: 已校验的 `cfg`、奇数 `gridSize`、方法名 `measureName`。
- Produces: `grid = uav.spectrum.upaIntegrationGrid(cfg,gridSize,measureName)`；`grid.nodes`为 `N×2`方向余弦，`grid.manifold`为 `M×N`，`grid.weights`为 `1×N`，并包含 `representedSolidAngle`、`measureName`和排除区掩码诊断。

- [ ] **Step 1: 写配置与网格失败测试**

在 `TestTask8UpaIncmDiagnosis.m` 中建立：

```matlab
function configurationContainsBothSnapshotCounts(testCase)
    cfg = config("smoke");
    d = cfg.literature.guLeshem2012.upaDiagnosis;
    testCase.verifyEqual(d.snapshotCounts, [8, 64]);
    testCase.verifyEqual(d.gridSizes, [21, 51, 101, 201]);
    testCase.verifyEqual(d.measureNames, ...
        ["direction_cosine", "solid_angle", "az_el_solid_angle"]);
end

function integrationGridsUsePositiveFiniteWeights(testCase)
    cfg = config("smoke");
    for name = cfg.literature.guLeshem2012.upaDiagnosis.measureNames
        grid = uav.spectrum.upaIntegrationGrid(cfg, 51, name);
        testCase.verifySize(grid.manifold, [cfg.derived.M, size(grid.nodes, 1)]);
        testCase.verifyTrue(all(isfinite(grid.weights)));
        testCase.verifyGreaterThan(min(grid.weights), 0);
        testCase.verifyFalse(any(grid.includedMask & grid.desiredMask));
    end
end
```

- [ ] **Step 2: 运行失败测试**

Run:

```text
matlab -batch "startup; r=runtests('tests/TestTask8UpaIncmDiagnosis.m'); assertSuccess(r)"
```

Expected: FAIL，因为 `upaDiagnosis`和 `upaIntegrationGrid`尚不存在。

- [ ] **Step 3: 添加唯一配置入口**

在 `config.m` 的 `cfg.literature.guLeshem2012` 中增加：

```matlab
'upaDiagnosis', struct( ...
    'snapshotCounts', [8, 64], ...
    'gridSizes', [21, 51, 101, 201], ...
    'measureNames', ["direction_cosine", "solid_angle", "az_el_solid_angle"], ...
    'selectedGridSize', 101, ...
    'selectedMeasureName', "solid_angle", ...
    'loadingRatio', 1e-2, ...
    'fitRegularization', 1e-4, ...
    'peakThresholdRatio', 0.15, ...
    'peakNmsRadiusBins', 2, ...
    'maximumCandidateCount', 8, ...
    'jitterSweepDeg', [0, 4, 8], ...
    'monteCarloTrials', literatureTrials, ...
    'periods', workload.sweepPeriods, ...
    'experimentIndex', 9)
```

在 `validateGuLeshem2012.m` 中逐项校验快拍数、奇数网格、唯一方法名、加载、正则化、候选数、抖动点、试次数、周期数和实验编号。

- [ ] **Step 4: 实现三类求积网格**

`upaIntegrationGrid`必须：

```matlab
function grid = upaIntegrationGrid(cfg, gridSize, measureName)
%UPAINTEGRATIONGRID 生成UPA二维诊断所需的节点、流形和求积权重。
```

- `direction_cosine`：在 `[-1,1]^2` 单元中心筛选可见圆盘，权重为 `du*dv`；
- `solid_angle`：使用同一单元中心节点，权重为 `du*dv/abs(uz)`；
- `az_el_solid_angle`：用方位角和高度角单元中心，权重为 `cos(el)*daz*del`，再转换为方向余弦；
- 三类方法都调用 `uav.array.steeringUPA`生成流形，并利用名义三维方向与配置排除半宽形成 `desiredMask`；
- 返回的 `includedMask`只表示排除期望区后的干扰积分节点。

- [ ] **Step 5: 运行Task 1测试并提交**

Run:

```text
matlab -batch "startup; r=runtests('tests/TestTask8UpaIncmDiagnosis.m'); assertSuccess(r)"
```

Commit:

```text
git add config.m src/+uav/+config/validateGuLeshem2012.m src/+uav/+spectrum/upaIntegrationGrid.m tests/TestTask8UpaIncmDiagnosis.m
git commit -m "实现UPA双快拍配置与积分网格"
```

---

### Task 2: 非负Hermitian协方差拟合核心

**Files:**
- Create: `src/+uav/+incm/fitNonnegativePowers.m`
- Modify: `tests/TestTask8UpaIncmDiagnosis.m`

**Interfaces:**
- Consumes: `Rx[M×M]`、`candidateManifold[M×Q]`、非负 `noisePower`、非负 `lambdaFit`。
- Produces: `[powers,diagnostics] = uav.incm.fitNonnegativePowers(Rx,candidateManifold,noisePower,lambdaFit)`；`powers[Q×1]>=0`，诊断包含拟合残差、字典秩和是否使用真值（固定为false）。

- [ ] **Step 1: 写两源功率恢复失败测试**

```matlab
function nonnegativeFitRecoversSyntheticSourcePowers(testCase)
    cfg = config("smoke");
    xi = [0.2, -0.4; -0.3, 0.1];
    A = uav.array.steeringUPA(xi, cfg.array.Nx, cfg.array.Ny, ...
        cfg.array.spacingOverLambda);
    expected = [20; 7];
    noise = 1;
    Rx = A * diag(expected) * A' + noise * eye(cfg.derived.M);
    [actual, diagnostics] = uav.incm.fitNonnegativePowers(Rx, A, noise, 0);
    testCase.verifyEqual(actual, expected, 'RelTol', 1e-9);
    testCase.verifyLessThan(diagnostics.relativeResidual, 1e-10);
    testCase.verifyFalse(diagnostics.usedTruth);
end
```

- [ ] **Step 2: 运行失败测试**

Expected: FAIL，函数不存在。

- [ ] **Step 3: 实现保Frobenius范数的Hermitian向量化和NNLS**

对Hermitian矩阵使用：实对角元素保持不变，严格下三角的实部和虚部分别乘 `sqrt(2)`。先构造 `residualCovariance=(Rx+Rx')/2-noisePower*I`，然后：

```matlab
augmentedDictionary = [dictionary; sqrt(lambdaFit) * eye(Q)];
augmentedTarget = [target; zeros(Q, 1)];
powers = lsqnonneg(augmentedDictionary, augmentedTarget);
```

拒绝非有限输入、负噪声、负正则化和空候选字典。

- [ ] **Step 4: 运行测试并提交**

Commit:

```text
git add src/+uav/+incm/fitNonnegativePowers.m tests/TestTask8UpaIncmDiagnosis.m
git commit -m "实现非负Hermitian协方差功率拟合"
```

---

### Task 3: `upa_covfit_incm`重构器与欠定快拍支持

**Files:**
- Create: `src/+uav/+incm/upaCovfit.m`
- Modify: `tests/TestTask8UpaIncmDiagnosis.m`

**Interfaces:**
- Consumes: 不含 `truth` 的 `obs`白名单和已校验 `cfg`。
- Produces: `[RinHat,diagnostics] = uav.incm.upaCovfit(obs,cfg)`；诊断包含加载量、噪声来源、期望候选、干扰候选、拟合功率、残差、矩阵条件数和 `usedTruth=false`。

- [ ] **Step 1: 写真值隔离与L=8失败测试**

构造零抖动64快拍观测并截取前8列：

```matlab
safe = uav.experiment.algorithmObservation(obs8);
testCase.verifyFalse(isfield(safe, 'truth'));
[RinHat, diagnostics] = uav.incm.upaCovfit(safe, cfg8);
testCase.verifyTrue(all(isfinite(RinHat), 'all'));
testCase.verifyLessThan(norm(RinHat-RinHat', 'fro'), 1e-10);
testCase.verifyGreaterThan(min(real(eig(RinHat))), 0);
testCase.verifyEqual(diagnostics.noiseSource, "configured-underdetermined");
testCase.verifyFalse(diagnostics.usedTruth);
```

- [ ] **Step 2: 运行失败测试**

Expected: FAIL，`upaCovfit`不存在。

- [ ] **Step 3: 实现候选提取、噪声估计和重构**

实现顺序固定为：

1. 把 `upaDiagnosis.loadingRatio`临时映射到 `cfg.spectrum.diagonalLoading`，调用 `loadedCapon`；
2. 在期望排除区内选择当前最强谱点，并始终加入名义期望导向矢量作为期望拟合字典；
3. 在排除区外用 `selectPeaks`提取最多 `min(maximumCandidateCount,M-1)`个干扰候选，无峰时回退到区外最强有限谱点；
4. `L<M`时使用 `cfg.signal.noisePower`；`L>=M`时使用较小 `ceil(M*noiseEigenvalueFraction)`个样本特征值的中位数；
5. 调用 `fitNonnegativePowers`同时拟合期望与干扰候选；
6. 丢弃期望候选功率，仅用干扰候选和 `noisePower*I`形成 `RinHat`；
7. 只做数值Hermitian化，不调用原有 `projectStructure`。

- [ ] **Step 4: 运行测试并提交**

Commit:

```text
git add src/+uav/+incm/upaCovfit.m tests/TestTask8UpaIncmDiagnosis.m
git commit -m "实现支持欠定快拍的UPA协方差拟合INCM"
```

---

### Task 4: 配对快拍构造与实验9方法执行器

**Files:**
- Create: `src/+uav/+scenario/truncatePeriodObservation.m`
- Create: `src/+uav/+experiment/runUpaDiagnosisMethods.m`
- Modify: `tests/TestTask8UpaIncmDiagnosis.m`

**Interfaces:**
- Produces: `obsShort = uav.scenario.truncatePeriodObservation(obsLong,L)`，重新计算 `Rx`并同步截取 `X`、`pilot`、真值干扰符号。
- Produces: `record = uav.experiment.runUpaDiagnosisMethods(obs,cfg,methodNames)`，方法固定允许 `oracle`、`esp`、`wco`、`gu_leshem2012_rec`和 `upa_covfit_incm`。

- [ ] **Step 1: 写前缀一致性和原始Rec不适用测试**

```matlab
function shortObservationIsExactPrefixAndOriginalRecIsUnavailable(testCase)
    [cfg64, obs64] = zeroJitterObservation(testCase, 64);
    obs8 = uav.scenario.truncatePeriodObservation(obs64, 8);
    testCase.verifyEqual(obs8.X, obs64.X(:, 1:8));
    testCase.verifyEqual(obs8.Rx, obs8.X*obs8.X'/8, 'AbsTol', 1e-13);
    cfg8 = cfg64; cfg8.timeline.L = 8; cfg8 = uav.config.validate(cfg8);
    record = uav.experiment.runUpaDiagnosisMethods(obs8, cfg8, ...
        ["oracle", "gu_leshem2012_rec", "upa_covfit_incm"]);
    testCase.verifyTrue(isnan(record.sinrLinear(2)));
    testCase.verifyEqual(record.status(2), "singular-not-applicable");
    testCase.verifyTrue(isfinite(record.sinrLinear(3)));
end
```

- [ ] **Step 2: 运行失败测试**

Expected: FAIL，两个新接口不存在。

- [ ] **Step 3: 实现观测截取和方法执行**

`truncatePeriodObservation`必须保留 `truth.trueIncm`等评价真值，但非Oracle方法进入 `algorithmObservation`白名单。执行器对 `L<M`的原始Rec不调用重构函数，保存NaN、状态和原因；其他方法共享同一 `obs.Rx`、名义导向和计时口径。

- [ ] **Step 4: 运行测试并提交**

Commit:

```text
git add src/+uav/+scenario/truncatePeriodObservation.m src/+uav/+experiment/runUpaDiagnosisMethods.m tests/TestTask8UpaIncmDiagnosis.m
git commit -m "实现UPA诊断的配对快拍与方法执行器"
```

---

### Task 5: 实验9编排、诊断指标和结果导出

**Files:**
- Create: `src/+uav/+experiment/runUpaIncmDiagnosis.m`
- Create: `experiments/exp09_upa_incm_diagnosis.m`
- Modify: `tests/TestTask8UpaIncmDiagnosis.m`

**Interfaces:**
- Produces: `result = exp09_upa_incm_diagnosis(cfg)`，包含 `raw.zeroJitter`、`raw.jitterSweep`、`summary`、方法状态、配置、元数据与六类结果文件。

- [ ] **Step 1: 写紧凑集成失败测试**

测试配置使用1次试验、1周期、网格 `[21,31]`、两种测度和 `L=[8,64]`，断言：

```matlab
result = exp09_upa_incm_diagnosis(cfg);
testCase.verifyEqual(result.snapshotCounts, [8,64]);
testCase.verifySize(result.raw.zeroJitter.sinrLinear, [1,2,2,2,5]);
testCase.verifySize(result.raw.jitterSweep.sinrLinear, [1,1,3,2,5]);
testCase.verifyTrue(all(result.raw.zeroJitter.snapshotPrefixMatched, 'all'));
testCase.verifyFalse(any(uav.experiment.entryPoints()== ...
    "exp09_upa_incm_diagnosis"));
```

并检查MAT、CSV、FIG、PDF、PNG和JSON全部存在。

- [ ] **Step 2: 运行失败测试**

Expected: FAIL，实验9入口不存在。

- [ ] **Step 3: 实现零抖动诊断数据流**

每个试次只生成一次64快拍观测，再截取8快拍。对每个网格和测度计算直接Rec积分诊断；`L=8`原始Rec保持NaN状态。保存迹比、最大特征值比、归一化特征值谱误差、条件数、子空间距离和真值干扰邻域评价指标。

- [ ] **Step 4: 实现抖动扫描和图表**

抖动扫描使用配置选定的网格/测度，同时运行标准MVDR版 `upa_covfit_incm`和同一INCM配合WCO权值的诊断曲线。图至少包含：

1. L=8与L=64的零抖动输出SINR；
2. 网格收敛；
3. 归一化特征值谱误差；
4. 0°/4°/8°抖动SINR；
5. 扇区违约率。

- [ ] **Step 5: 运行测试并提交**

Commit:

```text
git add src/+uav/+experiment/runUpaIncmDiagnosis.m experiments/exp09_upa_incm_diagnosis.m tests/TestTask8UpaIncmDiagnosis.m
git commit -m "实现UPA INCM诊断实验9"
```

---

### Task 6: Smoke、Preliminary性能闭环与参数锁定

**Files:**
- Modify only if diagnostics support it: `config.m`
- Modify: `docs/gu_leshem2012_rec_reproduction.md`
- Modify: `docs/experiment_protocol.md`
- Modify: `README.md`

**Interfaces:**
- Consumes: 实验9原始MAT结果。
- Produces: 固定配置下可复现的smoke/preliminary结果和中文诊断结论。

- [ ] **Step 1: 运行smoke反馈环**

```text
matlab -batch "startup; cfg=config('smoke'); cfg.plot.visible=false; cfg.plot.closeAfterExport=true; cfg.output.overwrite=true; r=exp09_upa_incm_diagnosis(cfg);"
```

逐项检验H1至H5的预测，不能只看SINR。

- [ ] **Step 2: 只依据诊断调整有物理含义的公共参数**

允许调整的只有加载比例、拟合正则化、峰阈值、NMS半径和候选数上限；同一组参数同时用于L=8与L=64。每次只改变一个参数，并保存评分表。禁止调整对照方法参数。

- [ ] **Step 3: 锁定参数并运行preliminary**

```text
matlab -batch "startup; cfg=config('preliminary'); cfg.plot.visible=false; cfg.plot.closeAfterExport=true; cfg.output.overwrite=true; r=exp09_upa_incm_diagnosis(cfg);"
```

验收零抖动L=8与L=64距Oracle均不超过2 dB，并检查51与101网格差异不超过0.5 dB。若失败，回到对应假设继续诊断，不得隐藏结果。

- [ ] **Step 4: 更新中文文档并提交**

文档必须报告实际试次、周期、快拍、网格、锁定参数、成功或失败的假设、精确SINR差距以及结论边界。

Commit:

```text
git add config.m README.md docs/experiment_protocol.md docs/gu_leshem2012_rec_reproduction.md
git commit -m "记录UPA INCM诊断与双快拍结果"
```

---

### Task 7: 最终回归、静态检查与独立审查

**Files:**
- Modify only for confirmed defects found by verification.

**Interfaces:**
- Produces: 可审查提交、完整测试证据和保留的结果文件。

- [ ] **Step 1: 运行专项测试**

```text
matlab -batch "startup; r=runtests('tests/TestTask8UpaIncmDiagnosis.m'); assertSuccess(r)"
```

- [ ] **Step 2: 运行全部测试并区分既有失败**

```text
matlab -batch "startup; r=runtests('tests','IncludeSubfolders',true); disp(r)"
```

预先存在的 `TestTask4TemporalRobust/temporalAblationsAreModuleLocalAndPilotDefaultsOff`必须单独报告；任何新增失败必须修复。

- [ ] **Step 3: 对全部新增或修改的MATLAB文件运行 `checkcode`**

要求新增文件0条静态检查消息；`git diff --check`无空白错误。

- [ ] **Step 4: 独立只读代码审查**

审查重点为真值隔离、L=8原始Rec状态、多测度Jacobian、协方差拟合向量化、配对随机性和结论边界。修复所有Critical与Important问题。

- [ ] **Step 5: 提交最终修复并保留分支**

```text
git status --short --branch
git log --oneline --decorate -5
```

不自动合并或推送；把分支、提交、测试结果和实验结果交给用户决定下一步。
