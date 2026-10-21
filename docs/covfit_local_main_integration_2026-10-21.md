# CovFit改名与局部细化版本接入main

## 本次范围

按用户要求，将研究分支`codex/lin-yue-upa-snapshot-study`的`1477a46`中
`refined_nnls0`接入main，供实验6、实验7使用。不合并整个研究分支，不移植
Le全分区细化、Lin/Yue论文核心或研究脚本，不修改开题报告。

|当前方法标识|含义|参数入口|
|---|---|---|
|`CovFit`|原`static_incm_mvdr`，原实现和默认参数不变|`cfg.incm.staticCovfit`|
|`CovFit_local`|本次移植的CovFit＋局部细化|`cfg.incm.localCovfit`|
|`static`|已有二维角域积分INCM，保持不变|原静态积分参数|

新标识区分大小写。旧字符串`static_incm_mvdr`作为兼容输入，在配置校验和
经典方法分发时规范为`CovFit`，新结果/图例不再输出旧名。若同时选旧名和
`CovFit`，按同义重复拒绝。历史MAT及历史结果文件不修改、不重新标记。
配置字段`staticCovfit`保留名称，避免破坏历史参数结构；它不是原积分`static`。

## 移植算法及不变项

新方法依次执行：

1. 以加载0.001的Capon谱建立候选和估计噪声；
2. 干扰区域采用八邻域局部极大值筛选，而不是最强相邻网格单元；
3. 在初始干扰峰附近作连续Capon峰细化，半径为2个方向余弦网格间隔；
4. 保留原期望拟合字典，联合进行非负协方差功率拟合，L2正则为0；
5. 删除期望组，只保留干扰项与噪声形成INCM；
6. 使用与原CovFit相同的固定名义期望导向和标准MVDR接口。

不使用真实干扰数、真实方向、真实姿态或真实INCM。允许的配置期望先验
仍用于固定分区。期望区观测峰只是拟合字典的一部分，不替换最终参考。

独立参数默认值：

|字段|默认值|说明|
|---|---:|---|
|`loadingRatio`|0.001|仅新方法的候选/细化前端加载|
|`fitRegularization`|0|最终NNLS正则，不改变旧CovFit的5000|
|`peakThresholdRatio`|0.15|相对候选阈值|
|`peakNmsRadiusBins`|0|候选筛选NMS半径|
|`maximumCandidateCount`|7|同时受M−1上限约束|
|`radiusBins`|2|局部搜索半径，不是角度2°|
|`kktTolerance`|1e−8|功率拟合返回点的数值检查阈值|
|`solver.maxIterations`|500|连续优化迭代上限|
|`solver.optimalityTolerance`|1e−9|最优性容差|
|`solver.constraintTolerance`|1e−6|优化器约束容差|
|`solver.stepTolerance`|1e−12|步长容差|
|`solver.feasibilityTolerance`|1e−8|候选接受时重新检查的可行性容差|

网格仍由`cfg.spectrum.gridSize`控制。main原生smoke默认21网格，而研究
smoke读取的是历史MAT中的51网格；不能只因都叫smoke就期待曲线绝对值一致。
一致性检查必须固定实际物理参数、网格、周期及随机子流。

新增或修改模块：

- `src/+uav/+incm/localCovfit.m`：独立新方法入口；候选/噪声接口复用旧函数，未使用的旧拟合结果不进入新INCM。
- `src/+uav/+spectrum/localPeakCandidates.m`：局部极大值候选筛选。
- `src/+uav/+spectrum/refineInterfererCandidates.m`：研究版数学步骤保留，参数路径改为main的`incm.localCovfit`，不依赖研究分支配置。
- `src/+uav/+incm/nonnegativeFitKkt.m`：独立数值检查，不裁剪或暗中修复功率。
- `runClassicMethods.m`：改名与新方法分发，均固定名义导向、标准MVDR。
- `runClassicSweep.m`：新名称、独立线形/标记与图例。
- `config.m`与`validate.m`：独立参数、可选组合、合法性与历史别名。
- 相关旧名称测试和已停用的历史调参入口同步规范名称；不重新启用积分调参协议。

KKT失败会明确报错，不静默回退旧CovFit，也不以裁剪负功率掩盖异常。研究
分支已观察到候选漏检、分区和联合拟合能量误分配造成的低性能样本；这次
移植不宣称解决这些边界，更不保证所有抖动强度都优于其他方法。

## 用户配置和运行

实验6/7默认七方法均包含两个CovFit；仍可自行改为子集：

```matlab
startup;
cfg = config("preliminary");
cfg.experiments.classicJitterMethodNames = ["oracle","CovFit","CovFit_local"];
cfg.experiments.classicSnapshotMethodNames = ["oracle","CovFit","CovFit_local"];
% 只调整新方法，旧CovFit不受影响。
cfg.incm.localCovfit.radiusBins = 2;
a = exp06_classic_method_jitter_sweep(cfg);
b = exp07_classic_method_snapshot_sweep(cfg);
```

常规完整实验仍输出`results/{raw,tables,figures,logs}`下的实验6/7原命名，
相同profile会覆盖同名结果。本次轻量核对输出到独立子目录，不覆盖用户的
正式preliminary结果。

## 小规模核对结果

仅运行`tests/TestCovfitLocalTransfer.m`，4项检查全部通过，没有运行完整
测试集、全量抖动/快拍扫描或WCO性能实验：

1. 规范名称、默认选择和历史别名/重复检查；
2. 非法局部半径被拒绝；
3. 从研究分支smoke归档复放L=8、32、1016各前2周期，两个CovFit共12条
   逐周期结果与研究值的线性SINR绝对差异**全部为0**；旧别名权值与新名一致。
   同时检查固定名义响应、标准MVDR与KKT；
4. 实验6/7各运行1试次×2周期×2扫描点，仅Oracle及两个CovFit，导出成功，
   全部方法共享相同快拍。

对照CSV：`results/tables/covfit_local_transfer_smoke/research_main_parity.csv`。
实验6/7的轻量FIG/PDF/PNG、MAT、CSV、日志分别保存在
`results/{figures,raw,tables,logs}/covfit_local_transfer_smoke/`。
这些轻量图只证明接入流程可运行，不替代研究分支20×20的性能结论。

核对命令：

```matlab
startup;
r = runtests("tests/TestCovfitLocalTransfer.m");
assertSuccess(r);
```

复放测试依赖本机研究分支保留的smoke MAT；在其他机器缺少该夹具时会明确
跳过该项，不能将跳过当作已通过一致性验证。
