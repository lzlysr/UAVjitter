# 实验7：ULA六方法快拍扫描smoke记录

## 配置及执行

- 分支：`codex/gu-leshem2012-rec-study`。结果生成时源码为未提交工作树，随后分职责归档提交；不改变main。
- 入口：`experiments/exp07_ula_classic_snapshot_sweep.m`；新增配置位于`cfg.experiments.ulaClassicSnapshot`。
- 阵列：1×9，沿y轴，间距半波长。ULA可观测角为theta=asin(u_y)，不使用UPA二维扫描。
- 扫描：MATLAB表达式`8:100:1024`，实际为[8,108,208,308,408,508,608,708,808,908,1008]，不包含1024。
- smoke计算量：3个独立试次、每点每次2周期；周期1ms。不同L点采用独立子流，同一点内六方法共享同一观测。
- 三轴姿态稳态RMS均为5°，rho=[0.95,0.95,0.98]；期望世界方向[0,-35]°，干扰方向[40,-25]°、[-55,-20]°，SNR=10dB，INR均为30dB。
- 保留动态姿态，不是零抖动或固定方向偏差实验。快拍数改变的是单周期样本数，不是时间长度。
- ULA一维积分与候选网格均为2001个中点节点，区间[-pi/2,pi/2]；期望排除区以名义可观测角为中心、半宽25°。
- WCO：用户确认采用文献平方半径比例0.3，即epsilon=sqrt(0.3×9)=1.6431676725。未使用导致不可行的原4°设计抖动分位数半径；不改变真实5°姿态分布，不修改WCO求解器。
- CovFit：复用已有`ulaCovfit`，不开启方向细化。参数取本分支`literature.guLeshem2012.upaDiagnosis`的锁定配置：加载0.12、L2正则5000、门限0.15、NMS半径0、候选上限7。它不是无正则网格CovFit研究变体。

## 方法名与适用性

顺序为Oracle、DL、ESP、WCO、Gu–Leshem原始Rec、CovFit。

本分支旧`runClassicMethods`中的`static_incm_mvdr`仍调用旧积分static。本次专项实验独立将该名称用于用户在main中的含义，即一维CovFit移植；未全局重写旧经典调度器，也不声称完全等同于main的二维候选流程。

原始Rec直接复用`guLeshem2012Rec`，采用未加载SCM的一维Capon积分，不扣底、不显式加白噪声项、不修复重构特征值。L=8小于9阵元，SCM必然秩亏，6个周期样本均标记为不适用，SINR保留NaN；其余快拍点成功。未使用伪逆冒充原始方法。

六方法最终使用同一个名义导向和标准MVDR接口，只有Oracle使用真实方向和真实INCM。原有MVDR接口自身的最小数值稳定加载规则保持不变。

## 结果与边界

本轮成功导出FIG/PDF/PNG、原始MAT、汇总CSV及JSON清单，均以`exp07_ula_classic_snapshot_sweep_smoke`命名，位于工程`results`的对应标准目录。全程约7.65秒。入口检查了有限结果与Oracle上界；不适用点单独记状态和原因，不伪装为数值零。

部分输出SINR中位数（dB）：

| L | Oracle | DL | ESP | WCO | 原始Rec | CovFit |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 19.40 | 3.97 | 4.31 | 3.13 | 不适用 | -1.46 |
| 108 | 19.45 | 8.34 | 17.24 | 10.30 | 17.47 | 17.90 |
| 508 | 19.40 | 12.05 | 18.37 | 16.04 | 17.63 | 18.66 |
| 1008 | 19.42 | 9.59 | 18.58 | 13.40 | 17.59 | 18.22 |

结果并非CovFit在所有快拍数下领先：8快拍时明显较差；1008快拍时ESP也高于CovFit。原始Rec在可逆SCM条件下接近17—18dB。smoke每点只有6个相关周期样本，不同扫描点的随机实现又不同，不能把曲线局部起伏唯一归因为L，也不能将本图推广为所有ULA几何下的方法排序。

## 运行命令

在本分支工作树中运行：

```matlab
startup;
r = exp07_ula_classic_snapshot_sweep(config("smoke"));
```

修改扫描点、阵元数、ULA网格、WCO平方半径比例，均在`config.m`中的`cfg.experiments.ulaClassicSnapshot`设置。其他场景参数仍读根配置。命令行批处理可增加`-singleCompThread`并在运行时关闭图窗，不影响数学定义。
