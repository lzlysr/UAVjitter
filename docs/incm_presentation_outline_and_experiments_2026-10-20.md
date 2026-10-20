# 三页组会汇报：前两页大纲与第三页实验运行说明

## 一、资料与结论边界

- 原理依据：`D:/workspace/hithesis-open/INCM原理深入思考-从SCM到INCM的根因.md`。
- 背景浓缩自：`D:/大组会/260917组会1.pptx`第2页，仅保留通信侧姿态扰动问题，不混入导航防欺骗。
- 原图：`D:/workspace/hithesis-open/examples/hitart/reports/figures/matlab原始图/32快拍20扰动范围下的鲁棒方法对比.fig`。
- 已读取原图图例，方法为理想上界、DL、ESP、WCO、INCM，横轴0—20°、21个点。该FIG的UserData为空，不能仅凭文件名和曲线恢复全部历史参数。新结果采用当前main完整配置，不宣称逐元素复现历史曲线。
- 当前只交付大纲与MATLAB实验支持，不生成PPTX，也不修改开题报告或上述原始FIG。

## 二、第1页：姿态抖动下的波束失配与INCM重构动机

### 页面布局

上方约三分之一介绍背景，下方用两组公式说明动机，页底保留一句结论。

### 背景文字

阵风、湍流和机体振动引起横滚、俯仰与偏航变化，改变阵列坐标系中的期望与干扰来向。名义期望响应与干扰零陷因此失配。缩短处理周期有利于局部方向稳定，却减少可用快拍，增加协方差估计误差。

目标：在姿态扰动与有限快拍下保持期望响应、抑制干扰泄漏，提高输出SINR。

### 核心公式

第一组定义观测统计量。这里省略周期下标，L为单周期快拍数。

$$
\widehat{\mathbf R}_x=\mathbf X\mathbf X^H/L,
\qquad
\mathbf R_x=P_s\mathbf a_s\mathbf a_s^H+\mathbf R_{i+n}.
$$

第二组在名义约束下展开目标。

$$
\min_{\mathbf w^H\bar{\mathbf a}_s=1}
\left[\mathbf w^H\mathbf R_{i+n}\mathbf w
+P_s|\mathbf w^H\mathbf a_s|^2\right].
$$

用颜色突出最后一项，旁注：方向失配时，降低真实期望响应也能降低目标代价。

页底结论：重构INCM，移除目标函数中的期望信号惩罚；INCM本身不自动校正期望方向。

### 口头解释与必要限定

- 上述目标分解使用准确总体协方差和源间不相关假设，不是有限快拍SCM的逐样本恒等分解。
- 方向准确时，约束使期望项成为常数，不产生上述额外自消诱因。
- 有限快拍仍会引入样本交叉项、病态性与噪声放大，不能把这些全部称为自消。
- 不在这一页塞入白化匹配系数和SINR闭式比值，避免符号过多。

背景文献可沿用旧PPT：Chen等，*Adaptive Hybrid Beamforming for UAV mmWave Communications Against Asymmetric Jitter*，IEEE TWC，2024，23(8):9432–9445，DOI `10.1109/TWC.2024.3362384`。

## 三、第2页：谱积分与候选功率拟合两条INCM路线

### 共同前置步骤

混合快拍、计算SCM、名义期望区域先验、加载Capon二维扫描、划分期望区与其补集。

这是流程层面的共同点，不代表当前两实现使用完全相同的加载矩阵：static读取`spectrum.diagonalLoading=0.01`，CovFit读取`incm.staticCovfit.loadingRatio=0.12`。当前默认对照保留各自历史参数。

### 左侧：现有static的二维谱积分

从区域划分后分叉：保留期望区外的稠密谱，扣除启发式谱底，对允许网格统一求积，加白噪声项，做结构处理。

$$
\widehat{\mathbf R}_{i+n}^{\rm static}
=\mathcal P\left[
\sum_{g\in\Omega_i}\omega_g\widetilde P_C(\boldsymbol\xi_g)
\mathbf a_g\mathbf a_g^H+\widehat\sigma_n^2\mathbf I
\right].
$$

其中ω为方向余弦单元面积，残余谱为扣底截断后的Capon谱，P表示现有Hermitian化、特征值下限及可配置收缩。无须在幻灯片上另展开所有子公式。

原理关注点：Capon峰值与峰体积不同；测度与背景谱会改变重构的相对功率和空间结构。数值积分收敛不等于真实INCM估计准确。

文献：[Luo等，*URGLQ: An Efficient Covariance Matrix Reconstruction Method for Robust Adaptive Beamforming*，IEEE TAES，2023，59(5):5634–5645，DOI `10.1109/TAES.2023.3263386`](https://doi.org/10.1109/TAES.2023.3263386)。正文第5636页式（10）—（12）明确Capon积分路线，后续用期望分量投影及Gauss–Legendre求积改善。**原文阵列为ULA，只作积分原理依据，不能称为普通UPA二维积分的直接验证，也不能把现有static称为URGLQ复现。** 2012年Gu–Leshem文献可放讲稿作为历史来源，不充当近年文献。

### 右侧：当前CovFit的候选功率拟合

区域划分后挑少量候选。期望区采用名义点及区内谱峰，干扰区执行阈值、峰值筛选和数量限制。各候选提前带有期望/干扰区域标签。

$$
\widehat{\mathbf p}
=\arg\min_{\mathbf p\ge0}
\left\|\widehat{\mathbf R}_x-\widehat\sigma_n^2\mathbf I
-\mathbf A_c\operatorname{diag}(\mathbf p)\mathbf A_c^H\right\|_F^2
+\lambda\|\mathbf p\|_2^2.
$$

图中用文字写最终合成：期望与干扰候选联合拟合，删除期望标签对应项，只累加干扰候选外积并加同一噪声项。Capon提供候选位置，不直接提供最终拟合功率，不乘求积面积。

原理关注点：候选漏检、误分组、离网误差、字典相关性与正则偏置。拟合SCM不代表消除了SCM误差。

文献：[Yue等，*Adaptive Beamforming for Cascaded Sparse Diversely Polarized Planar Array*，IEEE TVT，2023，72(12):15648–15664，DOI `10.1109/TVT.2023.3288611`](https://doi.org/10.1109/TVT.2023.3288611)。第15655页式（46）—（50）直接支持“联合拟合功率、仅保留干扰项并加噪声”。原文使用稀疏不同极化平面阵和闭式最小二乘，当前普通UPA、Capon候选和非负L2正则是工程变体，不冒称原样复现。两篇2023年论文按约近三年选用，不把2021文献伪写为近年。

### 底部汇合

同一名义参考导向、标准MVDR权值、真实环境评价输出SINR。两路均不加入期望导向校正和时序预测。汇报结论只写当前场景实测差异，不写UPA二维积分必然失效。

## 四、第3页的两张图与代码变化

图1：固定L=32，扫描三轴相同RMS为0:1:20°。

图2：固定三轴RMS为`cfg.attitude.stdDeg`，默认[5,5,5]°，扫描每周期快拍数。当前main配置为`8:100:512`，实际节点为[8,108,208,308,408,508]；下文保留的10×10历史预览使用[8,16,32,64]。横轴表示每个处理周期可用样本数，不是处理周期数，也不是累计训练长度。

默认两图均选择：Oracle、DL、ESP、WCO、原static、CovFit。分别由两个可编辑方法数组控制，并不写死。

修改范围：

- `config.m`：经典抖动扫描加入static，增加独立的`classicSnapshotMethodNames`。
- `src/+uav/+config/validate.m`：允许经典列表使用static，校验新增方法选择字段。
- `src/+uav/+experiment/runClassicMethods.m`：新加static分支，直接调用原重构器，再用标准MVDR；原static重构代码及主线鲁棒权值调用不变。
- `src/+uav/+experiment/runClassicSweep.m`：复用原实验6执行流程，统一两种扫描的配对观测、簇bootstrap统计及绘图导出。
- `runClassicJitterSweep.m`：保留原接口，转发至共用流程。
- `runClassicSnapshotSweep.m`及`experiments/exp07_classic_method_snapshot_sweep.m`：新增快拍扫描入口，不替换原实验3，也不加入run_all。
- 两图图例区分积分INCM（static）和CovFit-INCM，Oracle保持黑色虚线，不固定纵轴裁切范围。

两方法最终权值接口相同，但它们自身的谱加载、噪声估计与正则设置并不相同；完整算法优劣不能全部归因于某一个步骤。本次不搜索参数，也不刻意降低其他方法性能。

## 五、怎样运行及修改参数

在`D:/workspace/UAVjitter`中运行：

```matlab
startup;
cfg = config("preliminary");
jitterResult = exp06_classic_method_jitter_sweep(cfg);
snapshotResult = exp07_classic_method_snapshot_sweep(cfg);
```

正式参数在`config.m`集中修改：

- 图1快拍数：`cfg.timeline.L`，本次保持32。
- 图1横轴：`cfg.experiments.jitterSweepDeg`。
- 图2固定抖动：`cfg.attitude.stdDeg`，默认[5,5,5]；该字段不会控制图1的扫描点。
- 图2横轴：`cfg.experiments.snapshotSweep`。
- 图1方法：`cfg.experiments.classicJitterMethodNames`。
- 图2方法：`cfg.experiments.classicSnapshotMethodNames`。
- 试次：`cfg.monteCarlo.nTrials`；两图各点周期数分别为`periodCounts.jitterSweep`、`periodCounts.snapshotSweep`。preliminary当前为50次、每点50周期。
- 原static参数：`cfg.spectrum`、`cfg.incm`；CovFit专用参数：`cfg.incm.staticCovfit`；DL/ESP/WCO参数位置不变。

需要快速核对执行流程时将profile改为smoke。smoke只有3试次、每点2周期，其曲线不能当作可靠的性能结论。

MATLAB交互式单实验默认显示并保留图窗。命令行`matlab -batch`进程结束后图窗不能继续存在，可通过保存的FIG重新打开。

## 六、输出位置

默认工程`results`下：

- `figures/exp06_classic_method_jitter_sweep_preliminary.{fig,pdf,png}`。
- `figures/exp07_classic_method_snapshot_sweep_preliminary.{fig,pdf,png}`。
- `raw/`下对应同名MAT，保留完整配置、原始逐周期指标和随机子流。
- `tables/`下对应CSV，包含中位SINR、95%中位数置信区间及5%低分位。
- `logs/`下对应JSON发布清单。

smoke文件以`_smoke`结尾。相同profile重跑会按现有`overwrite=true`规则替换工程内同名结果，不修改开题报告目录中的原始图。程序运行时打印每个扫描点完成情况，避免长时间没有进度提示。

## 七、本次10×10预览的复现方式

用户确认先采用10次试验、每点10周期的预览规模。仅在运行时配置副本覆盖计算量，`config.m`中的preliminary默认50×50保持不变。51网格、32快拍抖动扫描、固定三轴5°的快拍扫描，以及所有物理/算法参数均保持默认。完整50×50运行已停止，未将中途结果发布为完整实验。

```matlab
startup;
cfg = config("preliminary");
cfg.monteCarlo.nTrials = 10;
cfg.experiments.periodCounts.jitterSweep = 10;
cfg.experiments.periodCounts.snapshotSweep = 10;
cfg.output.rawDir = fullfile(cfg.output.rawDir, 'preview_10x10');
cfg.output.figureDir = fullfile(cfg.output.figureDir, 'preview_10x10');
cfg.output.tableDir = fullfile(cfg.output.tableDir, 'preview_10x10');
cfg.output.logDir = fullfile(cfg.output.logDir, 'preview_10x10');
cfg = uav.config.validate(cfg);
jitterResult = exp06_classic_method_jitter_sweep(cfg);
snapshotResult = exp07_classic_method_snapshot_sweep(cfg);
```

输出分别位于`results/{raw,figures,tables,logs}/preview_10x10`。虽然文件名沿用`_preliminary`以表示51网格等参数来源，但实际试次和周期数在MAT配置和日志元数据中明确保存为10×10，不得当作50×50结果引用。命令行运行增加`-singleCompThread`只控制计算线程数，不改变算法或扫描参数。

### 当前运行状态与WCO阻断记录

2026-10-20，本轮10×10抖动扫描21个点已经完成并发布。快拍扫描在首点L=2因WCO数值优化失败而停止，**尚未发布完整的10×10快拍扫描图**；不能把此前smoke图冒充本轮结果。

固定子流复现定位为实验7、扫描点1、试次3、周期10。SCM秩为2，最小特征值约−6.32×10^−13，谱范数约9397.14；WCO半径约1.84432，小于参考导向范数3。原fmincon返回exitflag=−3。仅将目标协方差除以正尺度的只读诊断虽然返回exitflag=2，但权值范数约9.19×10^6，按原SCM计算二次型约−12.82，不能据此认定有效求解。

本轮未修改WCO模型、半径、正则或数值求解器，未更换种子，也未删除低快拍点。失败复现输入暂存于被Git忽略的`data/cache/preview_wco_failure.mat`。需先确定数值求解修复或失败点记录策略，再继续交付第二张图。

随后用户明确要求快拍扫描从8开始，因此将默认扫描点改为[8,16,32,64]，不实施WCO求解器修复。仍使用10×10运行时规模和原物理/算法参数，重新运行实验7。随机子流继续按扫描点序号分配，新列表中的L=8对应点1，不与旧列表中点3的L=8宣称逐元素配对。

重新运行已完成，四个快拍点均成功，所有方法的原始SINR均为有限值，运行约32.35秒。FIG/PDF/PNG和MAT/CSV/JSON均已发布至上述`preview_10x10`目录。此次只重跑实验7，已完成的实验6图不变。完整默认试次与周期数仍为50×50，WCO求解器未修改。

| 快拍数L | Oracle | DL | ESP | WCO | static | CovFit |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 19.39 | 5.63 | 5.59 | 5.73 | 4.91 | 13.41 |
| 16 | 19.29 | 8.31 | 8.41 | 8.61 | 4.90 | 14.45 |
| 32 | 19.37 | 11.49 | 11.72 | 11.67 | 5.31 | 15.21 |
| 64 | 19.37 | 13.04 | 13.99 | 13.65 | 3.94 | 15.40 |

表中为输出SINR中位数（dB），仅作本次10×10预览结果，不代表所有几何或抖动条件下的方法排序。每个点内方法配对，但不同快拍点采用不同子流；不能将曲线局部波动唯一归因于快拍数。

### 后续main配置调整与历史结果区分

本次提交前，用户将`snapshotSweep`调整为`8:100:512`，并将随机子流的`experimentStride`从10^8调整为10^9，同时补充了试次、周期及bootstrap参数的解释注释。这些改动保留入库，不更改上述历史预览数据。新配置重跑时的快拍节点与随机子流编号均已不同，因此不能宣称本提交下直接运行会逐元素复现上述表格；复核旧图应读取其原始MAT中保存的完整配置与种子。
