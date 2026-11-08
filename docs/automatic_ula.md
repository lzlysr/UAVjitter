# main的ULA自动降维与两条INCM路线

## 实现范围与参考来源

- `Nx>1 && Ny==1`自动识别为`ULA`，其余既有尺寸保留历史路径。本次不扩展y轴ULA或单阵元算法。
- `Nx>1 && Ny>1`保持原UPA网格、参数与重构定义。main默认仍为3×3，不自动改成9×1。
- 参考`codex/gu-leshem2012-rec-study`的`d666399`：`steeringULA`、`ulaStatic`、`ulaCovfit`、`guLeshem2012Rec`及ULA专项exp7。
- 该分支专项exp7实际上使用y轴ULA和独立参数。本次仅借鉴一维等角中点积分与拟合思路，改为x轴投影，不复制其2001点网格、实验方向、全分区细化或期望方向估计。
- 姿态三轴过程、世界方向、信号功率、快拍数和随机流继续读取main。没有增加期望导向校正，最终参考仍是固定名义导向。

## 模型与积分约定

三维链路仍为`姿态 -> 旋转 -> 阵列坐标方向 -> UPA/ULA导向 -> 快拍`。
ULA观测的是`u_x`，等效角为`theta=asin(u_x)`，不是世界方位角。
扫描节点为`theta_j=-pi/2+(j-1/2)*pi/G`，`G=cfg.spectrum.gridSize`。
方向节点统一表示为`[sin(theta_j),0]`。第二坐标置零只是等效表示，不是估计到了真实的另一方向维度。

ULA求积权重为`pi/G`，即角度弧度长度`dtheta`。UPA继续使用方向余弦面积`dxi_x*dxi_y`。
不得把这两种测度误称为相同的角域密度；改变测度需要独立实验，不能暗中归一化功率。
兼容旧代码的`cellArea`在ULA中表示单元角长度；新诊断显式记录`dimension`、`coordinate`、`measure`、`step`和`weights`。

期望排除区在ULA中为`abs(theta-theta_nominal)<=ulaDesiredExclusionHalfWidthDeg`（计算时转换为弧度）。
UPA继续读取`desiredExclusionHalfWidthDeg`。统一取参函数为`uav.spectrum.desiredExclusionHalfWidth`，
实际生效宽度记录于`cfg.derived.desiredExclusionHalfWidthDeg`、算法诊断和结果元数据。
它排除整个相同x轴投影的响应等价类，不能用不可观测的y坐标重新分类期望能量。
谱为行向量；所有候选节点强制为N×2，拟合流形为M×N，避免MATLAB行向量索引导致意外转置或隐式广播。

## 方法定义保持分离

| 配置名称 | 积分/拟合定义 | 最终经典实验权值 |
| --- | --- | --- |
| `incm_integral_origin` | 未加载SCM的Capon谱直接求积；不扣底、不显式加噪声、不做PSD或收缩修复 | 固定名义导向的标准MVDR |
| `incm_integral_processed` | 加载Capon、扣谱底、均匀求积、显式噪声及原结构处理 | 同一标准MVDR |
| `incm_covfit` | 网格候选、期望/干扰联合非负功率拟合、只合成干扰项及噪声 | 同一标准MVDR |
| `incm_covfit_refined` | 局部干扰峰、连续一维局部细化、联合非负拟合；不细化期望字典 | 同一标准MVDR |

历史输入名`gu_leshem2012_rec`、`static`、`CovFit`、`CovFit_local`依次映射为表中四个新名称；
`static_incm_mvdr`映射为`incm_covfit`。映射保留原顺序，但拒绝同义重复。
配置、调度、MAT/CSV、诊断及图例统一输出新名称；不重写已存在的历史结果文件。
内部函数名`guLeshem2012Rec`、`upaCovfit`、`localCovfit`及参数结构`staticCovfit`、`localCovfit`保留，以减少历史配置断裂。

原式Rec核心直接来自研究分支，接口为`guLeshem2012Rec(Rx,A,weights)`。
ULA的原式积分使用`dtheta`；UPA若选择该名称，是原式在现有二维网格/测度下的迁移，不声称复现论文原始阵列场景。
L<M或SCM数值奇异时，原始Rec报专用异常；经典调度器只捕获这一种已知不适用情形，记录NaN及原因，其余异常原样抛出。
不使用伪逆或偷偷加载把它变成另一个算法。统计仅使用有效样本，CSV记录有效数量与不适用数量/原因。

局部细化的`radiusBins`在UPA中仍对应方向余弦网格步长；ULA中对应角网格步长`pi/G`。
因此不同维度下同样的bin数不意味着同样的物理搜索半径；一维细化诊断额外保存`radiusRad`。
`temporal`也支持一维谱平移与一维自适应求积；这里只保证功能适配，不宣称时序先验因此得到验证。

## ULA专用鲁棒半径

实测main原设计分布的95%误差半径在9×1 ULA下不可行。经用户授权，ULA采用研究分支已有的固定设置：

`epsilonS=sqrt(cfg.beamforming.calibration.ulaSquaredRadiusScale*M)`，默认比例0.3。

M=9时半径约1.64317，小于名义导向范数3。必须满足`0<scale<1`。
该固定设计只依赖配置与阵元数，不依赖当前试次、扫描点或性能。UPA继续原来的离线分位数标定。
主线鲁棒权值共用该标定入口；经典exp6/7中的各INCM方法仍用标准MVDR，不获得额外WCO半径收益。
经典ULA实验未选择WCO时跳过无关的半径标定。

## 文件职责

- `src/+uav/+config/validate.m`：识别ULA、建立一维角网格、派生元数据及名义几何诊断。
- `src/+uav/+spectrum/admissibleDirections.m`：一维可观测期望排除区。
- `src/+uav/+spectrum/loadedCapon.m`：复用同一矩阵求解，输出一维谱及测度。
- `src/+uav/+spectrum/selectPeaks.m`、`localPeakCandidates.m`：单行网格自然只比较左右邻居；候选节点固定N×2。
- `src/+uav/+spectrum/refineUlaCandidates.m`：只优化干扰的一维角度，保持候选局部边界及排除区。
- `src/+uav/+spectrum/ulaAdaptiveQuadrature.m`：一维细分及粗块聚合，保持允许网格单元的角长度权重总和。
- `src/+uav/+spectrum/align2D.m`：ULA将行平移范围固定为零，只搜索列平移；保留原无证据回退。
- `src/+uav/+incm/reconstruct.m`、`guLeshem2012Rec.m`：两种积分定义分开调度。
- `src/+uav/+incm/upaCovfit.m`：沿用拟合核心，仅将ULA名义候选转成唯一等效表示。
- `src/+uav/+experiment/runClassicMethods.m`：所有经典方法共享观测；原式Rec不适用时不阻断其他方法。
- `src/+uav/+experiment/collectClassicMethods.m`、`runClassicSweep.m`：逐周期状态、缺失原因、有效样本统计及新方法图例。
- `src/+uav/+experiment/finalize.m`、`metadata.m`：ULA文件名后缀和阵列/积分测度元数据。
- `tests/TestAutomaticUla.m`：降维、独立求积、欠采样、候选、时序、权重守恒及exp6/exp7轻量导出。

## 运行与输出

```matlab
startup
cfg=config("smoke");
cfg.array.Nx=9; cfg.array.Ny=1;
names=["oracle","dl_mvdr","esp","wco","incm_integral_origin","incm_integral_processed","incm_covfit","incm_covfit_refined"];
cfg.experiments.classicJitterMethodNames=names;
cfg.experiments.classicSnapshotMethodNames=names;
cfg=uav.config.validate(cfg);
disp(cfg.derived.ulaGeometry)
r6=exp06_classic_method_jitter_sweep(cfg);
r7=exp07_classic_method_snapshot_sweep(cfg);
```

快拍扫描改`cfg.experiments.snapshotSweep`，抖动扫描改`jitterSweepDeg`，网格点数改`cfg.spectrum.gridSize`。
改完后实验入口会重新validate；手工调用底层算法时必须自行validate更新缓存。
main的smoke/preliminary/paper网格点数原样保留，并未自动换为研究分支2001点，所以这里只是公式/流程移植，不代表原论文完整数值复现。
增加网格密度时应研究求积收敛，不以单个网格结果断定路线失效。

标准结果位置仍为`results/raw`、`results/figures`、`results/tables`和`results/logs`。
ULA文件例如`exp07_classic_method_snapshot_sweep_ULA_9x1_smoke.fig`；UPA文件名保持原样。
针对性测试的轻量结果单独输出至上述目录的`automatic_ula_smoke`子目录，不覆盖用户已有正式实验。

## 几何与结论边界

main当前期望世界方向[0,-35]度、干扰[40,-25]及[-55,-20]度，在名义姿态下的x轴等效角约为55、43.97、32.61度。
最近名义角间隔约11.03度。ULA现在使用独立的5度半宽，所以名义干扰不再落入期望区；UPA的25度不变。
若把ULA宽度设回25度，会重新出现两个干扰均被排除的问题。`derived.ulaGeometry`记录名义重叠供用户判断，算法不能用它识别真实源身份。
这是显式、可记录的分区参数调整，不是算法数学形式的改进；也不是依据SINR筛选出的最优宽度。
宽度过小会漏掉抖动后的期望能量，过大又会纳入干扰。5度半宽不等于5度RMS的95%覆盖范围。
动态场景中还可能发生期望和干扰相对固定分区的交叉，无法靠一个固定宽度保证所有周期都正确分组。
ULA默认5度与UPA默认25度代表不同的排除先验；直接跨阵列比较时，不能把性能差异全部归因于阵列维数，需显式报告并做排除宽度敏感性分析。
因此建议先在零抖动下做路线/网格诊断，再在固定几何与配对随机流下扫描排除半宽及抖动强度，分别报告泄漏、漏干扰和输出SINR。
缺少ULA专用字段的历史MAT配置回退到原公共半宽（通常25度），不静默改为新默认5度；要使用新方案，请显式补充ULA字段。
相同阵元数的9×1与3×3仍有不同孔径和源响应相关性，不能把SINR差异全归因为积分维数。
应依次检查零抖动理论协方差、网格收敛、已知候选方向的Oracle诊断以及有限快拍/动态姿态，区分数值积分误差、谱功率偏差和方向先验失配。

## 首轮自动降维验证记录（2026-11-08，独立半宽与改名前）

- `TestAutomaticUla`新增9项测试，通过一维流形/权重、等价响应排除、独立原式求积、8快拍原始Rec不适用、64快拍有限性、空干扰、零抖动、一维时序、固定半径及轻量exp6/exp7导出检查。
- 扩展回归共67项：65项通过、2项失败。几何场景11/11通过，基线与静态INCM 21/21通过；既有时序测试21/22通过，CovFit移植测试3/4通过。
- 两个失败分别为`temporalAblationsAreModuleLocalAndPilotDefaultsOff`和`renamedMethodsAndLegacyAlias`。在隔离临时目录中解包未修改main的HEAD `a6a887a`重跑这两项，同样失败；不是本次新增回归。
- 前者在main收缩参数为0时仍要求结构消融造成大于1e-12的变化，实际只有约6.5459e-14；后者仍假定默认组合同时包含普通CovFit，但用户已经从默认组合移除了它。本次不改用户参数，也不改这两个历史测试。
- 历史研究数据的逐周期UPA数值复核仍通过：`incm_covfit_refined`的线性SINR逐项一致，普通CovFit差异约1e-12量级，满足原相对容差。
- 轻量导出仅每点1试次×1周期，抖动点[0,2]度、快拍点[8,64]，用于接口检查，不具有算法排序或学术统计结论的证据强度。没有运行preliminary全量扫描。

## 改名与独立排除区的验证范围

`tests/TestIncmNamesAndSectors.m`覆盖新旧名称在UPA/ULA下的数值一致性、统一重构接口、
同义重复拒绝、ULA的5/25度分区变化、UPA不受ULA宽度影响、非法参数拒绝，以及旧MAT宽度回退。
`TestAutomaticUla`的接口/导出测试使用新名称和新默认ULA半宽。旧方法名的兼容测试有意保留旧字面值。
随改名更新了历史方法名单断言，移除了“默认组合必须包含普通拟合”的过时假设；没有把普通拟合重新强加到默认组合。

本轮针对性验证44项全部通过：新名称与分区7项、ULA自动降维9项、经典积分/拟合扫描3项、CovFit历史移植复核4项、基线与静态INCM21项。
补充经典调度/图例/默认名单/单实验导出4项也全部通过，合计48项；这不是整个仓库全量测试通过的声明。
只运行了轻量导出与数值回归，没有运行新的preliminary性能扫描，也没有把分区宽度5度宣称为最优参数。
