# 细化 CovFit 的时序增强

## 方法入口和信息流

`incm_covfit_refined_temporal` 是实验 6、7、8 可选的 INCM 重构方法。
每周期先得到原 `incm_covfit_refined` 的连续方向细化结果和联合非负拟合结果。
两种方法一起运行时复用该计算，运行时间仍分别计入完整的单周期计算成本。

新方法的状态只包含上一周期已估计的有效干扰方向余弦和拟合功率。
每个蒙特卡洛试次重新初始化。首周期、无历史或没有通过关联时，
直接用当前周期的细化 CovFit。关联成功并完成非负拟合时，直接输出时序 INCM。
ULA 沿用单周期细化 CovFit，本次验证针对 UPA。最终权值采用与基线相同的名义期望
导向和标准 MVDR。时序重构使用当前 SCM、阵列配置及上一周期的干扰状态，
真实姿态、来向、INCM 及干扰数用于场景生成或指标计算。

## 每周期方向更新和跨周期关联

每个周期内同一路干扰仍对应一个点源。姿态改变使这个点在阵列坐标系中移动。
功率连续性采用已经匹配的历史系数，而导向矢量始终由当前周期的方向重算。
这样不会直接平均两个指向不同方向的 INCM，也不会把过去的点硬放到本周期的同一网格上。

设 $\mathbf U_k$ 是从当前 SCM 中抽出的显著特征子空间。
对 UPA 的 $x$、$y$ 两轴分别选择相邻的等形子阵，构造

$$
\widehat{\boldsymbol\Psi}_{r,k}
=(\mathbf J_{r,1}\mathbf U_k)^\dagger
 \mathbf J_{r,2}\mathbf U_k,\qquad r\in\{x,y\}.
$$

两个方向的旋转算子共享信号模态特征基。第 $j$ 个模态在第 $r$ 轴的复相位
为 $2\pi(d_r/\lambda)\xi_{r,j,k}$，由此估计两个方向余弦。然后在候选流形
$\mathbf A(\boldsymbol\xi)$ 上局部减小加权子空间的投影残差

$$
\left\|
(\mathbf I-\mathbf A\mathbf A^\dagger)\mathbf U_k\mathbf W_k^{1/2}
\right\|_F^2,\qquad
[\mathbf W_k^{1/2}]_{jj}
=\frac{\ell_j-\widehat\sigma_{n,k}^2}{\sqrt{\ell_j}},
$$

其中 $\ell_j$ 为对应 SCM 特征值。候选数由相对于估计噪声的特征值门限和子阵维数决定。
若当前方向候选与上一周期的某个干扰在方向余弦距离和候选功率比上均相容，
就按距离做一对一关联。关联成功的干扰用当前方向替代旧位置。
与期望拟合字典太接近的候选不纳入干扰组。较强且尚未被关联结果覆盖的当前
CovFit 候选可以进入干扰组，使新出现的干扰仍能被纳入拟合。

## 联合拟合与 INCM

以 $\mathcal M_k$ 表示关联成功的当前干扰集合，$\pi_k(j)$ 表示第 $j$ 路当前干扰
对应的上一周期索引。当前期望候选和干扰候选仍联合估计功率。目标函数为

$$
\min_{\mathbf q\ge0}
\left\|\widehat{\mathbf R}_{x,k}
-\sum_\ell q_\ell\mathbf a_{\ell,k}\mathbf a_{\ell,k}^H
-\widehat\sigma_{n,k}^2\mathbf I\right\|_F^2
+\lambda_0\|\mathbf q\|_2^2
+M^2\lambda_t\sum_{j\in\mathcal M_k}
\big(q_j-\widehat q_{\pi_k(j),k-1}\big)^2 .
$$

$\lambda_0$ 是细化 CovFit 共用的正则参数。$\lambda_t$ 来自
`cfg.incm.temporalCovfit.powerContinuityWeight`。这里采用二次连续性项，
可以直接给原非负最小二乘增广，即新增一行
$\sqrt{M^2\lambda_t}(q_j-\widehat q_{\pi_k(j),k-1})$。单个阵元幅度为 1 时，
导向矢量外积的 Frobenius 范数平方为 $M^2$，用它给时序项定标。
未匹配的干扰和期望候选不受历史功率约束。实现对应原研究计划中的
“关联后约束干扰功率”，将候选方案的绝对差换成了平方差。

拟合后只将干扰组组装为 INCM

$$
\widehat{\mathbf R}_{i+n,k}
=\sum_{j\in\mathcal I_k}\widehat q_{j,k}\,
\mathbf a_{j,k}\mathbf a_{j,k}^H
+\widehat\sigma_{n,k}^2\mathbf I.
$$

每个候选仍是单点外积，系数是点源功率。名义期望导向误差仍可能限制大抖动下
的输出 SINR，尤其应同时观察干扰泄漏、期望响应和 SINR 第 5 百分位。


## 验证口径

实验 6 扫描 $0{:}1{:}20^\circ$，每个非零点要求新方法中位数严格高于基线。
实验 7 扫描 $L=32{:}32{:}1024$，要求每个点均不低于基线。实验 8 扫描
输入 SNR $-10{:}5{:}30$ dB，逐点报告与基线的差值，单列退化点。
三组实验均用 3×3 UPA 和 `preliminary` 的 50 次试验、每次 50 个周期。
实验 7、8 的三轴 RMS 均为 5°。实验 6、8 的每周期快拍数为 32。
基线、时序方法和理想上界共享同一批快拍。逐点配对置信区间以完整试次为
bootstrap 重采样单位，每个子样本分别计算两种方法的中位 SINR，再取其 dB 差。

## 完整 preliminary 结果

本轮采用固定基种子 20260827、51 点每轴谱网格和同一组算法参数完成三种扫描。
每个扫描点和每种方法都有 2500 个有效周期样本。下表的增益为时序版中位数减
原细化版中位数，配对 95% 区间为逐点区间。

| 实验 | 中位 SINR 高于基线的点数 | 中位 SINR 增益范围 | 配对区间下限大于 0 的点数 |
|---|---:|---:|---:|
| 6，1°–20° 抖动 | 20/20 | +0.0353 至 +0.3838 dB | 20/20 |
| 7，32–1024 快拍 | 32/32 | +0.0026 至 +0.0408 dB | 32/32 |
| 8，−10 至 30 dB 输入 SNR | 7/9 | −0.0299 至 +0.0518 dB | 7/9 |

实验 6 的 0° 对照点提高 0.0393 dB。在 20° 点，输出 SINR 中位数从
16.318 dB 提高至 16.702 dB。

实验 8 的 −5 dB 和 0 dB 输入点分别下降 0.0299 dB 和 0.0198 dB，
其配对 95% 区间分别为 [−0.0406, −0.0200] dB 和 [−0.0483, −0.0004] dB。
这两个点的完整观测和置信区间与其他扫描点一同保留。

当前方法的结果位于 `results/diagnostics/temporal_covfit_no_pilot_preliminary/results`，
包含逐周期 MAT、FIG/PDF/PNG、CSV 和 JSON 清单。此前
`results/diagnostics/temporal_covfit_preliminary` 中的结果对应旧的导频残差选择版本。
[实验 6 CSV](../results/diagnostics/temporal_covfit_no_pilot_preliminary/results/tables/exp06_classic_method_jitter_sweep_preliminary.csv)、
[实验 7 CSV](../results/diagnostics/temporal_covfit_no_pilot_preliminary/results/tables/exp07_classic_method_snapshot_sweep_preliminary.csv)、
[实验 8 CSV](../results/diagnostics/temporal_covfit_no_pilot_preliminary/results/tables/exp08_classic_method_snr_sweep_preliminary.csv)。

## 时序贡献的消融

在实验 6 的 1°、10° 和 20° 三点，用相同的试次、周期和快拍进行四组对照。
仅当前周期组使用当前基线候选作为关联参考，采用相同的子空间细化并将历史功率项置零。
历史关联组改用上一周期的候选，仍关闭功率项。完整时序组同时启用历史关联和
功率连续性。这能分开比较历史与新增单周期处理带来的收益。

| 各轴 RMS | 原细化 CovFit | 仅当前周期 | 历史关联 | 完整时序 | 完整时序较仅当前周期的配对 95% 区间 |
|---|---:|---:|---:|---:|---:|
| 1° | 19.072 dB | 19.124 dB | 19.127 dB | 19.127 dB | [约 0, +0.0090] dB |
| 10° | 18.427 dB | 18.502 dB | 18.515 dB | 18.516 dB | [+0.0020, +0.0331] dB |
| 20° | 16.318 dB | 16.416 dB | 16.686 dB | 16.702 dB | [+0.066, +0.774] dB |

20° 点的完整时序相对仅当前周期组提高约 0.286 dB。只开历史关联已达到
16.686 dB，说明这一点的额外收益主要来自方向关联。1° 点的额外时序收益
约为 0.003 dB，整体收益主要来自新增的单周期子空间处理。20° 的第 5 百分位
仍为 −36.18 dB，大抖动下的深衰落仍需单独评价。

[消融数据](../results/diagnostics/temporal_covfit_no_pilot_preliminary/exp06_temporal_ablation.csv)和
[消融脚本](../results/diagnostics/temporal_covfit_no_pilot_preliminary/run_temporal_ablation.m)。


## 全域 Capon 与子空间前端的对照

实验 6、7、8 可以选择 `incm_covfit_refined_temporal_capon`。它从可见圆盘的加载
Capon 谱选择局部峰，再按连续 Capon 目标细化。期望排除区内外的峰均可被提名。
细化只受可见圆盘和粗网格峰附近的搜索范围限制。候选的干扰身份在此之后判定。

子空间版 `incm_covfit_refined_temporal` 继续用两轴子阵的平移不变性提取候选。
两条链路共用候选关联功率的双边最小二乘公式、距离和功率比匹配、期望候选的
近邻排除、联合功率拟合与标准 MVDR。它们各自维护上周期状态。Capon 谱峰高度
用于前端筛选，最终 INCM 中的系数来自联合协方差拟合。

当前周期的候选位置由当前数据估计，历史方向只用于匹配和决定哪些点进入干扰
字典。将预测方向或其协方差写入连续细化目标，作为后续方向先验研究路线。

完整对照采用相同的 3×3 UPA 场景、`preliminary` 的 50 次试验 × 50 个周期、
相同的快拍和随机子流。方法列表包含 Oracle、原细化 CovFit、子空间时序版及
全域 Capon 时序版。导出的 `caponVsSubspaceMedianDeltaDb` 为全域 Capon 版
中位数减子空间版中位数，置信区间按完整试次进行配对 bootstrap。

~~~matlab
startup
cfg = config("preliminary");
methods = ["oracle", "incm_covfit_refined", "incm_covfit_refined_temporal", ...
    "incm_covfit_refined_temporal_capon"];
cfg.experiments.classicJitterMethodNames = methods;
cfg.experiments.classicSnapshotMethodNames = methods;
cfg.experiments.classicSnrMethodNames = methods;
cfg.plot.closeAfterExport = true;
cfg.output = struct('rootDir', ...
    fullfile(pwd,'results','diagnostics','temporal_covfit_full_capon_comparison'), ...
    'overwrite', false);
r6 = exp06_classic_method_jitter_sweep(cfg);
r7 = exp07_classic_method_snapshot_sweep(cfg);
r8 = exp08_classic_method_snr_sweep(cfg);
~~~


### 完整配对结果

| 实验 | 子空间版中位 SINR 较高的点数 | Capon 减子空间的中位 SINR 差值范围 |
|---|---:|---:|
| 实验 6，1°–20° 抖动 | 20/20 | −0.3213 至 −0.0353 dB |
| 实验 7，32–1024 快拍 | 32/32 | −0.0397 至 −0.0026 dB |
| 实验 8，−10 至 30 dB 输入 SNR | 7/9 | −0.0511 至 +0.0303 dB |

0° 对照点另计，全域 Capon 版低 0.0393 dB。实验 8 中全域 Capon 版在输入
SNR 为 −5 dB 和 0 dB 时分别高 0.0303 dB 和 0.0206 dB。两处优势的配对 95%
区间分别为 [+0.0205, +0.0411] dB 和 [+0.0016, +0.0484] dB。上述 61 个点的
逐点配对区间均未跨零，其中包含只有几千分之一 dB 的差异。子空间版的中位数
与原无导频对照一致。结果适用于上述两套前端和当前 3×3 UPA 配置。

[实验 6 CSV](../results/diagnostics/temporal_covfit_full_capon_comparison/results/tables/exp06_classic_method_jitter_sweep_preliminary.csv)、
[实验 7 CSV](../results/diagnostics/temporal_covfit_full_capon_comparison/results/tables/exp07_classic_method_snapshot_sweep_preliminary.csv)、
[实验 8 CSV](../results/diagnostics/temporal_covfit_full_capon_comparison/results/tables/exp08_classic_method_snr_sweep_preliminary.csv)。
逐周期 MAT 及 FIG/PDF/PNG 与这些表保存在同一次运行的结果目录中。

### 20° 抖动点诊断

在实验 6 的 20° 点，原细化 CovFit、全域 Capon 时序版和子空间时序版的中位
SINR 分别为 16.318、16.381 和 16.702 dB。Capon 相对子空间的配对差值为
−0.321 dB，95% 区间为 [−0.816, −0.118] dB。

同一组观测另做候选覆盖诊断。每次试验从第 2 周期开始统计，共 2450 个周期。
两路真实干扰各自至少有一个候选在方向余弦距离 0.04 内时，记为覆盖。
真值只参与事后评价。关联前，Capon 与子空间的覆盖率分别为 96.98% 和 99.84%。
关联、筛选及功率拟合后，只计正功率干扰点，两版覆盖率分别为 70.78% 和 72.65%。
两版的输出干扰泄漏中位数分别为 −20.28 和 −21.07 dB。
这些数据支持子空间版在该点更容易把真实干扰对应的方向纳入模型。
这是对两条完整候选前端的比较，两者的候选检测和细化方式均有差异。


## 期望区内子空间、区外 Capon 的混合对照

`incm_covfit_refined_temporal_hybrid` 复用两种现有方向前端。按**细化后的方向**
使用原有期望区划分规则，保留区内的子空间候选和区外的 Capon 候选。两侧候选
落入同一去重邻域时，优先保留区内的子空间点。合并后在同一候选流形上重新估计
关联所需的功率，避免使用不同候选字典下的功率作为同一组匹配的依据。

三种时序版共用历史关联规则、功率连续性系数、期望拟合字典和标准 MVDR。
各版按试次独立保存历史状态。分区负责筛选候选，区内点仍须通过历史关联后才能
作为干扰纳入 INCM。历史方向只用于关联。将方向先验纳入细化目标仍是后续研究路线。

混合版在 3×3 UPA 的完整 `preliminary` 设置下与两种纯前端配对比较。
沿用上一节复现命令，把 `methods` 加上 `"incm_covfit_refined_temporal_hybrid"`，
并将结果根目录改为 `results/diagnostics/temporal_covfit_hybrid_comparison`。
CSV 的 `hybridVsSubspaceMedianDeltaDb` 和 `hybridVsCaponMedianDeltaDb`
分别记录混合版减子空间版、混合版减 Capon 版的中位 SINR 差。


### 完整混合对照结果

| 实验 | 混合版减纯子空间版 | 混合版减纯 Capon 版 |
|---|---:|---:|
| 实验 6，1°–20° 抖动 | −0.2999 至 −0.0353 dB，20 个点均较低 | −0.0027 至 +0.0214 dB |
| 实验 7，32–1024 快拍 | −0.0397 至 −0.0026 dB，32 个点均较低 | 32 个点的中位数均相同 |
| 实验 8，−10 至 30 dB 输入 SNR | −0.0512 至 +0.0303 dB，仅 −5 dB 和 0 dB 输入点较高 | 最大绝对差 0.0006 dB |

差值采用输出 SINR 中位数。实验 6 的 0° 点另计，混合版与纯 Capon 的中位数相同，
比纯子空间低 0.0393 dB。三组实验中，混合版相较纯 Capon 的逐点配对 95% 区间
均包含零或接触零。两条纯前端曲线与上一轮重跑结果一致。

[实验 6 CSV](../results/diagnostics/temporal_covfit_hybrid_comparison/results/tables/exp06_classic_method_jitter_sweep_preliminary.csv)、
[实验 7 CSV](../results/diagnostics/temporal_covfit_hybrid_comparison/results/tables/exp07_classic_method_snapshot_sweep_preliminary.csv)、
[实验 8 CSV](../results/diagnostics/temporal_covfit_hybrid_comparison/results/tables/exp08_classic_method_snr_sweep_preliminary.csv)。

### 20° 抖动点的分区诊断

排除每次试验的首周期后，共有 4900 个真实干扰来向实例，其中 665 个在期望区内，
4235 个在区外。单路干扰在候选集中有一个方向余弦距离不超过 0.04 的候选即记为
覆盖。真值及其区内外标签仅用于事后评价。算法按估计候选的方向划区。

| 诊断量 | 全子空间 | 全域 Capon | 分区混合 |
|---|---:|---:|---:|
| 关联前区内干扰覆盖率 | 98.95% | 98.20% | 98.80% |
| 关联前区外干扰覆盖率 | 99.98% | 98.44% | 98.39% |
| 拟合后两路均有正功率候选的周期比例 | 72.65% | 70.78% | 70.82% |
| 输出干扰泄漏中位数 | −21.07 dB | −20.28 dB | −20.23 dB |

在这个场景，86.4% 的真实干扰来向实例在区外。混合版提高小部分区内来向的
候选覆盖率，区外覆盖率仍接近全域 Capon。20° 的输出中位 SINR 分别为
16.702、16.381、16.402 dB。混合版相较子空间版的配对差值为 −0.300 dB，
95% 区间为 [−0.796, −0.117] dB。相比全域 Capon，混合版在该点的点估计
提高 0.021 dB，95% 区间下限约为 0。
