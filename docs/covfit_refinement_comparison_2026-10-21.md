# 两种CovFit细化的移植比较与32快拍参数扫描

## 固定实验合同

本轮在`codex/lin-yue-upa-snapshot-study`继续；main、Gu–Leshem研究分支和开题报告均不修改。使用`karpathy-guidelines`的最小必要改动原则，新增独立方法，不覆盖原CovFit或论文功率核心的定义。

所有非Oracle重构最终固定相同名义期望导向，接同一标准MVDR。算法不读取真实方向、真实干扰数、真实姿态和真实INCM。所有真值只在Oracle或事后指标评价层使用。

## 三个需要分别命名的版本

1. `refined_nnls0`：本次局部细化。51网格、低加载0.001、八邻域局部峰、局部半径2个方向余弦网格间隔，只细化干扰候选，期望拟合字典保持不变，NNLS正则为0。
2. `le_refined_covfit`：从Gu–Leshem研究分支`d666399`原样移植的完整细化版本。201网格、加载0.001、最强网格单元候选、全先验分区约束搜索、内部期望拟合候选与干扰候选均细化，细化后按1e−5方向余弦容差去重，NNLS正则为0。保留名义期望候选，最终期望参考不变。本轮排除区25°，不移植额外的7°窄区优势。
3. `le_shared_candidates`：机制对照，不是旧完整算法。把Le细化器接到本次相同的51网格局部峰；只细化干扰，保留相同期望字典和候选列数，不去重。用于控制输入，分析局部限制与全分区优化的差别。

移植文件为`upaCovfitLe.m`与`refineCaponNodesLe.m`，仅改函数名和配置路径，并增加观测真值字段拒绝检查；数学流程、目标、约束、去重和默认求解设置不变。已在相同3组输入（L=8、32、1016）上分别运行两个分支，三组重构矩阵相对误差均为0，拟合功率也一致。完整版本不等于Gu–Leshem原始Rec，仍是工程CovFit改进。

## 实验A：重跑快拍数图

`runLinYueSnapshotStudy`保留原24方法，追加以上两个Le版本；沿用原图归档20试次×20周期、0°/5°两条件、L=[8,24,40,72,136,264,520,1016]及原随机子流。WCO继续只引用同观测5°归档，不重新调参或改基线。

原先结果在`lin_yue_2026-10-21`中保留；本轮输出标签改为`lin_yue_refinement_2026-10-21`，不会覆盖上一轮MAT、图和审计。

## 实验B：32快拍单因素参数扫描

物理场景仍沿用原图归档：3×3、半波长、SNR10 dB、INR30/30 dB、三轴5°RMS，其余几何与姿态相关性不变。固定L=32；preliminary为20试次×20周期，smoke为3×2。

新实验使用实验7的未使用扫描点子流编号65；全部参数点与全部方法共用同一批观测，不为某个参数改变随机实现。四个参数分别扫描，其他项回到默认值，不做笛卡尔联合调优：

|参数轴|扫描值|作用范围|
|---|---|---|
|局部半径，网格间隔数|0.25、0.5、1、2、4、8|只影响本次局部细化及共享它的Lin/Yue核心；Le全分区方法为不受影响的参照|
|最大迭代次数|5、10、20、50、100、200|同时修改本次与Le细化器；未成功退出时按各自原规则保留初始节点|
|最优性容差|1e−4、1e−6、1e−8、1e−10|同时修改两种细化器；其他容差不变|
|谱加载比例|1e−4、1e−3、1e−2、0.12|修改两种细化前端；原CovFit基线仍保持原0.12配置|

加载属于候选/细化前端参数，不是求解器收敛容差；扫描它会改变谱和初始候选，不能把其效应解释成仅由局部优化产生。201与51网格、期望拟合候选处理、去重和各自其他容差默认不同，完整方法比较不能单独归因于搜索范围。

比较方法可在`config.m`的`linYueStudy.refinementSweep.methodNames`选择。默认包括Oracle、原CovFit、本次细化、Le完整细化、Le同候选细化、Yue-LS核心和Lin点功率核心。后两者仍是共同前端的核心移植，不是完整论文复现。

## 数值和统计处理

32快拍实验对每次非负拟合检查KKT残差，阈值1e−8；异常返回点标记失败，不裁剪功率、不伪造有效点。无约束LS重构非PSD也保留失败。局部优化未收敛的候选按原规则回退到初始节点，同时保存细化接受比例；这不等于整个波束形成失败。

所有方法保存逐周期SINR、期望增益、干扰泄漏、INCM误差、KKT、候选数量、细化接受比例、最终数值加载、状态与错误信息。SINR在线性域取中位数后转dB，95%区间按试次成簇bootstrap。含失败方法的曲线是成功子集统计，须同时查看CSV失败数。

`compareCovfitRefinements`同时报告Le相对本次局部细化的配对中位数比、簇bootstrap区间、逐周期差异5%/95%分位和最大绝对差异。0.1 dB只作事先明确的“数值接近”描述阈值，不凭不显著就宣称所有场景统计等效。参数扫描仅用于敏感性分析，不从中挑最好参数替换默认条件。

## 配置与运行

配置集中在`config.m`的`cfg.experiments.linYueStudy`：旧方法的`leGridSize/leCovfit/leRefinement`，本次的`localRefinement/refineCandidateRadiusBins`，新实验的`refinementSweep`。

```matlab
startup;
a = runLinYueSnapshotStudy("preliminary");
plotLinYueSnapshotStudy(a.artifacts.mat);
compareCovfitRefinements(a.artifacts.mat);
b = runCovfitRefinementSweep("preliminary");
compareCovfitRefinements(b.artifacts.mat);
```

输出统一在本工作树`results/{raw,tables,figures,logs}/lin_yue_refinement_2026-10-21/`。
快拍主图名为`lin_yue_preliminary_comparison`，32快拍参数图名为
`covfit_refinement_parameter_sweep_preliminary`，均导出FIG/PDF/PNG。
原始MAT包含实际配置、软件/代码版本、随机子流和逐周期状态；配对表另以`_refinement_paired.csv`结尾。

## 正式结果

正式运行后在此记录全部预定参数轴的结果、失败数和解释边界，不按smoke结果删除不利点。
