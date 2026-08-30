# 无人机姿态抖动下鲁棒 INCM 重构仿真系统

本工程是面向开题阶段初步研究工作的 MATLAB 仿真系统，用于研究无人机三轴姿态抖动如何改变机载矩形均匀平面阵列（UPA）的接收方向和阵列流形，并比较不同波束形成与干扰加噪声协方差矩阵（INCM）重构方法的输出性能。

当前场景为：无人机接收来自服务地面基站的 5G NR FR1 蜂窝下行期望信号，同时接收若干正常邻区地面基站产生的同频干扰。所有基站在世界坐标系中的方向固定，方向的时变性只由无人机阵列坐标系随姿态转动产生。

首版采用窄带、远场、单径等效基带模型，主要评价指标为波束形成器输出 SINR。工程暂不包含完整 5G NR 波形、多径宽带信道、基站运动、姿态传感器输入、学习型方向预测和最终预测—波束形成级联。

## 1. 软件环境

- MATLAB R2024a 或更高版本；
- Optimization Toolbox：用于求解球形导向矢量失配约束下的鲁棒权值；
- MATLAB Unit Test：用于运行自动化测试。

## 2. 快速开始

运行某个实验：
```matlab
cd('D:\workspace\UAVjitter')
startup

cfg = config("smoke");
cfg.plot.visible = true;

result = exp01_model_validation(cfg);
```

单独调用某个实验时，`cfg.plot.closeAfterExport` 默认为 `false`，因此可见
图窗在结果导出后继续保留。`run_all(...)` 属于批处理入口，会自动把该
开关设为 `true`，每组实验导出完成后立即关闭图窗，避免累计占用资源。

在 MATLAB 命令窗口中运行：

```matlab
cd('D:\workspace\UAVjitter')
startup

% 运行全部自动化测试
r = runtests('tests', 'IncludeSubfolders', true);
assertSuccess(r)

% 运行测试并执行五组小规模冒烟实验
run_all('smoke')
```

也可以从 PowerShell 直接运行：

```powershell
matlab -batch "cd('D:\workspace\UAVjitter'); startup; r=runtests('tests','IncludeSubfolders',true); assertSuccess(r)"
matlab -batch "cd('D:\workspace\UAVjitter'); startup; run_all('smoke')"
```

`config(profile)` 支持三种配置规模：

| profile | 用途 | 蒙特卡洛规模 |
|---|---|---:|
| `"smoke"` | 快速检查全部代码、实验和结果导出链路 | 3 次配对试验 |
| `"preliminary"` | 生成可用于初步分析的结果 | 100 次配对试验 |
| `"paper"` | 生成论文级统计结果 | 500 次配对试验 |

所有正式数值参数、扫描范围、随机种子、绘图设置和输出路径都集中在 `config.m`。实验脚本和算法文件不应自行写入场景数值。

## 3. MATLAB 目录名前为什么有“+”

`src/+uav` 以及其下的 `+array`、`+scenario` 等目录使用了 MATLAB 的包目录语法。目录名前的加号不是普通文件名装饰，而是 MATLAB 用来声明命名空间的固定规则。

例如文件：

```text
src/+uav/+array/steeringUPA.m
```

在程序中应调用为：

```matlab
a = uav.array.steeringUPA(xi, Nx, Ny, spacingOverLambda);
```

调用时不写加号。加号只出现在磁盘目录名中，其作用包括：

1. 把函数按职责组织到独立命名空间；
2. 避免不同模块出现同名函数时相互覆盖；
3. 使调用位置明确显示函数属于阵列、场景、INCM还是实验模块；
4. 避免把 `src` 下几十个子目录全部单独加入 MATLAB 路径。

因此，`+uav/+array` 不能简单重命名成 `uav/array`。删除加号后，`uav.array.steeringUPA` 这种调用方式将失效。

## 4. 系统数据流

单次仿真周期的数据流为：

```text
config(profile)
    ↓ 参数校验、单位转换和派生量计算
uav.config.validate
    ↓
uav.scenario.generate
    ↓ 生成三轴姿态、旋转矩阵和动态阵列流形
uav.scenario.generatePeriod
    ↓ 生成当前周期快拍 X、样本协方差 Rx 和隔离真值 truth
uav.experiment.runMethods / runVariants
    ↓
频谱估计 → INCM 重构 → 波束形成权值 → 指标计算
    ↓
统计汇总、绘图和结果文件导出
```

真值集中保存在 `obs.truth` 中，只能用于 Oracle 方法和性能指标计算。静态及时序 INCM 等非 Oracle 算法必须先经过 `uav.experiment.algorithmObservation` 生成白名单观测，算法只能获得 `X`、`Rx`、训练序列、周期索引和名义导向矢量。

## 5. 工程目录

```text
D:\workspace\UAVjitter\
├─ config.m                 唯一正式参数入口
├─ startup.m                将 src、experiments 和 tests 加入路径
├─ run_all.m                运行测试和五组主实验
├─ experiments\             五组主实验及三个独立研究实验入口
├─ src\+uav\               按 MATLAB 包组织的模型与算法源码
├─ tests\                   七个阶段的自动化测试
├─ docs\                    模型契约、变量字典和实验协议
├─ data\cache\              可重新生成的缓存
└─ results\                 MAT、CSV、FIG、PDF、PNG 和 JSON 结果
```

## 6. 根目录文件

### `config.m`

系统唯一参数入口。它负责设置：

- 阵列规模和阵元间距；
- 仿真周期、每周期快拍数；
- 姿态抖动稳态标准差和一阶高斯—马尔可夫相关系数；
- 期望信号与邻区干扰的世界坐标方向；
- SNR、INR和噪声功率；
- Capon谱、谱对齐、峰值选择和求积参数；
- INCM稀疏正则化与结构投影参数；
- 鲁棒权值和离线误差半径标定参数；
- 蒙特卡洛次数、扫描点、随机子流方案；
- 绘图和结果输出设置。

其中 `cfg.plot.fontSize` 统一控制所有实验图中的横轴标题、纵轴标题、
坐标轴刻度和图例字号，默认值为 `12`。

### `startup.m`

只负责加入运行所需路径，不生成参数、不运行实验，也不修改全局随机数状态。

### `run_all.m`

统一运行入口。它先执行全部测试，再按固定顺序调用五个实验，并返回各实验的结果文件路径。

## 7. `experiments`：主实验与独立研究实验入口

该目录中的文件都是很薄的用户入口。它们不包含数值参数，只负责接收配置或 profile，然后调用 `uav.experiment` 包中的正式实现。实验1至实验5由 `run_all` 批量执行；实验6至实验8是针对特定研究问题的独立实验，必须由用户显式调用。

### `exp01_model_validation.m`

建模与机理验证。输出三轴姿态、方向余弦和导向矢量失配随时间的变化，并检查零抖动场景下动态量是否保持不变。

### `exp02_jitter_sweep.m`

姿态抖动强度扫描。在每个扫描点的相同随机实现上比较所选方法的输出
SINR、期望方向增益和干扰泄漏，统计中位数、低分位数及配对置信区间；
同时计算与方法无关的期望信号归一化导向矢量失配，并以三轴共同抖动
RMS为横轴输出中位数及5%/95%试次分位范围。

实验2的方法组合由 `cfg.experiments.exp02MethodNames` 单独控制，可从
`oracle`、`frozen`、`scm`、`dl`、`static` 和 `temporal` 中选择任意
非空且不重复的子集，并保持配置顺序。例如：

```matlab
cfg = config("preliminary");
cfg.experiments.exp02MethodNames = ...
    ["oracle", "frozen", "scm", "static"];
cfg = uav.config.validate(cfg);
result = exp02_jitter_sweep(cfg);
```

该字段只影响实验2。实验3、实验4仍使用固定的
`cfg.experiments.methodNames`，独立经典方法实验仍使用
`cfg.experiments.classicJitterMethodNames`。

### `exp03_snapshot_sweep.m`

有限快拍扫描。改变每周期快拍数，比较输出 SINR、INCM归一化 Frobenius 误差、最小特征值、条件数和运行时间。

### `exp04_dynamic_stability.m`

动态稳定性实验。保留逐试次、逐周期结果，绘制 SINR 时序、经验CDF、期望增益、干扰泄漏及运行时间，并计算低于阈值的概率。

### `exp05_ablation.m`

时序方法消融实验。在完全相同的场景与快拍上依次关闭时间连续性、谱对齐、自适应求积和结构增强，并可开启训练序列投影。程序会验证每个变体只改变一个指定开关。

### `exp06_classic_method_jitter_sweep.m`

经典方法姿态抖动退化实验。方法组合由 `cfg.experiments.classicJitterMethodNames` 控制，用于研究经典波束形成器的输出SINR随三轴抖动RMS变化的趋势，不加入 `run_all`。

### `exp07_gu_leshem2012_reproduction.m`

Gu与Leshem（2012）论文例1的一维ULA忠实复现。实验分别扫描输入SNR和快拍数，固定论文给出的10阵元半波长ULA、期望方向5°、两个干扰方向−50°/−20°、每个干扰30 dB以及期望扇区[0°,10°]。该入口只比较本工程已有且原论文确实比较过的方法：理想先验MVDR、ESP、WCO以及 `gu_leshem2012_rec`。

### `exp08_gu_leshem2012_upa_transfer.m`

Gu–Leshem Rec从一维ULA直接迁移到当前矩形UPA的诊断实验。它只把一维角度积分改成二维方向余弦面积求积，不加入对角加载、噪声底座、结构投影、动态扇区、谱对齐或时序稀疏增强。实验同时输出SINR、INCM误差、主特征值比例、归一化特征值谱误差、条件数、干扰子空间距离和扇区违约率，用于分辨“原公式本身是否可复现”和“直接搬到UPA后为何退化”这两个问题。

## 8. `src/+uav`：模型与算法模块

### 8.1 `+config`：配置校验与派生量

- `validate.m`：检查配置字段、参数范围和维度；把角度转换到模型使用的弧度语义；计算阵列规模 `M`、干扰数 `J`、信号功率、世界方向、名义阵列流形以及二维Capon网格流形缓存；同时规范化输出目录。
- `validateGuLeshem2012.m`：独立校验文献复现与UPA迁移参数，限制方法集合为论文和当前实现的交集，并阻止未经说明的奇异样本协方差回退、偶数二维网格或缺失Rec方法。

### 8.2 `+geometry`：坐标系与姿态旋转

- `rotationZYX.m`：按照 `Rz(yaw)Ry(pitch)Rx(roll)` 生成阵列/机体系到世界坐标系的三维旋转矩阵，即满足 `v_world=R*v_array`，并检查输入维度。
- `worldToArray.m`：利用旋转矩阵的转置，把世界坐标系中的固定单位方向转换到当前阵列坐标系。
- `directionCoordinates.m`：从三维阵列坐标方向提取二维方向余弦，并计算用于展示和诊断的方位角、俯仰角。

### 8.3 `+array`：矩形UPA阵列流形

- `steeringUPA.m`：根据二维方向余弦、阵列行列数和波长归一化间距生成导向矢量；阵元顺序固定为 `vx⊗vy`，支持 `Nx≠Ny` 以及一次输入多个方向。
- `steeringULA.m`：按照论文的一维ULA角度定义生成导向矢量，采用 `exp(j2πdm sinθ)`；该函数只服务于文献忠实复现，避免把当前UPA的二维方向余弦约定误套到原论文上。

### 8.4 `+scenario`：动态场景与快拍生成

- `generate.m`：按照三轴独立一阶高斯—马尔可夫模型生成整段姿态轨迹，然后逐周期计算旋转矩阵、阵列坐标方向、方向余弦和真实导向矢量；不预生成大型快拍数组。
- `generatePeriod.m`：只生成指定周期的期望QPSK训练序列、独立干扰QPSK序列和复高斯白噪声，形成 `X`、`Rx`；真实INCM和真实阵列流形集中放入 `truth` 字段。
- `generateGuLeshem2012Trial.m`：生成论文例1的一次块快拍。期望、干扰和噪声均为独立复高斯过程；每次蒙特卡洛试验抽取固定的±4°实际DOA偏差，并构造期望扇区补集上的梯形求积节点和权重。

### 8.5 `+spectrum`：二维空间谱处理

- `loadedCapon.m`：对样本协方差进行尺度化对角加载，复用配置阶段缓存的二维可见域阵列流形，通过一次批量线性求解计算Capon功率谱。
- `admissibleDirections.m`：把“位于方向余弦可见圆盘内”和“位于期望信号排除区外”合并为连续方向可接受条件。
- `align2D.m`：在受限二维平移范围内计算当前谱与上一周期谱的归一化互相关；证据不足或候选歧义时退化为恒等对齐。
- `selectPeaks.m`：在增强谱上执行阈值筛选和非极大值抑制；峰数只受观测阈值和 `M-1` 上限限制，不使用真实干扰数。
- `adaptiveQuadrature.m`：峰值邻域采用细网格，其他干扰区域采用粗网格，并为方向余弦单元分配面积权重。
- `mergeQuadratureNodes.m`：合并细网格与粗网格中的重复或近重复节点，累计相应求积权重，避免同一区域重复积分。
- `guLeshem2012UpaQuadrature.m`：为UPA迁移实验生成期望排除区外的二维方向余弦节点及面积权重。它是“最直接二维化”的诊断模块，不声称等价于论文的一维角度积分。

### 8.6 `+incm`：INCM重构

- `reconstruct.m`：静态和时序INCM的统一入口。静态方法使用当前周期Capon谱和均匀求积；时序方法依次进行谱对齐、非负时序稀疏增强、峰值选择、自适应求积和结构投影。该文件还负责噪声功率估计、期望区排除和可选训练序列投影。
- `temporalProx.m`：求解非负逐元素时序稀疏目标。利用目标函数的分段凸结构枚举精确候选点，不依赖CVX。
- `projectStructure.m`：对原始INCM依次执行Hermitian化、半正定投影、特征值下限和可配置收缩，保证输出有限、Hermitian且半正定。
- `guLeshem2012Rec.m`：严格按论文式(6)计算Capon功率，再按式(7)对期望扇区补集积分。该函数不添加对角加载、伪逆、显式噪声单位阵、谱底座扣除、PSD修复或收缩；样本协方差不可逆时明确报错，而不是静默改变算法。

### 8.7 `+beamforming`：波束形成权值

- `solve.m`：统一求解MVDR、对角加载MVDR和球形导向失配约束下的鲁棒权值。鲁棒问题转换为 `2M` 维实变量后交给 Optimization Toolbox；当 `epsilonS=0` 时严格退化为标准MVDR。
- `calibrateEpsilon.m`：在固定设计抖动分布上离线生成导向误差样本，取配置分位数作为鲁棒误差半径。该标定与当前实验试次真值隔离，每个实验配置只执行一次。

### 8.8 `+metrics`：统一指标

- `evaluate.m`：使用真实期望导向矢量、真实INCM和源功率计算线性域输出SINR、期望阵列增益与干扰泄漏；dB转换在统计展示阶段完成，避免错误地直接平均dB值。

### 8.9 `+experiment`：实验编排与结果管理

#### 观测与方法执行

- `algorithmObservation.m`：从完整观测中构造非Oracle白名单，明确剥离 `truth` 及任何额外诊断字段。
- `runMethods.m`：在同一个周期快拍上依次运行Oracle、冻结、SCM、DL、静态INCM和时序INCM六种方法，并计算统一指标。
- `runVariants.m`：运行六个时序算法消融变体，并深度检查每个变体只改变指定配置字段。
- `collectMethods.m`：按“实验点—试次—周期”生成确定性子流，收集六种方法的逐试次逐周期原始指标。
- `collectVariants.m`：采用同样的配对随机机制收集消融实验原始指标。

#### 五个实验正式实现

- `runExp01.m`：建模验证的完整计算、绘图和导出逻辑。
- `runExp02.m`：抖动强度扫描的完整计算、统计和导出逻辑。
- `runExp03.m`：快拍数扫描的完整计算、统计和导出逻辑。
- `runExp04.m`：动态稳定性实验的完整计算、统计和导出逻辑。
- `runExp05.m`：模块消融的完整计算、配置诊断和导出逻辑。
- `entryPoints.m`：返回五个实验入口的固定顺序，供 `run_all` 和集成测试使用。

此外，`runClassicJitterSweep.m` 与 `experiments/exp06_classic_method_jitter_sweep.m`
构成独立的经典方法退化实验。该实验不加入 `run_all`，只在用户明确调用
时执行，并根据 `cfg.experiments.classicJitterMethodNames` 选择方法组合。

文献研究分为三个实现层：`runGuLeshem2012Methods.m`在同一观测上执行
Oracle、ESP、WCO和Rec；`runGuLeshem2012Reproduction.m`组织原论文的一维
SNR/快拍扫描；`runGuLeshem2012UpaTransfer.m`组织当前UPA场景的二维直接
迁移和扇区、网格诊断。后两者均为独立实验，不改变五组主实验入口。

#### 配置、随机性和统计汇总

- `prepare.m`：重新校验实验配置，并在实验开始时只标定一次 `epsilonS`。
- `substream.m`：根据实验、扫描点、试次、周期和用途分配互不重叠的确定性随机子流。
- `summarizeLinear.m`：对二维线性域样本计算中位数、低分位数和共享重采样置信区间。
- `summarizeTrials.m`：把同一试次内多个周期视为相关簇，以“试次”为独立bootstrap单位，避免把周期错误当成独立样本。
- `summarizePoints.m`：逐扫描点调用试次簇统计，用于抖动和快拍扫描。
- `summarizeFinite.m`：对可能含 `NaN` 的有限指标计算中位数和配置分位数。
- `summarizeFinitePoints.m`：逐扫描点汇总INCM误差、特征值、条件数和运行时间等有限指标。

#### 元数据、绘图和结果发布

- `metadata.m`：记录Git提交号、MATLAB版本、工具箱版本、随机种子方案、配置规模和生成时间。
- `finalize.m`：把配置与元数据加入结果并调用导出器；导出失败时始终关闭图窗，导出成功后根据 `cfg.plot.closeAfterExport` 决定是否关闭。单实验默认保留，`run_all` 批处理自动关闭。
- `writeArtifacts.m`：先在目标目录中生成完整临时结果包，再统一发布MAT、CSV、FIG、PDF、PNG和JSON。
- `publishStagedBundle.m`：实现带提交点的结果发布事务；提交前失败时恢复旧文件，提交后备份清理失败时保留完整新结果和未删除备份。

### 8.10 `+plot`：绘图辅助

- `newFigure.m`：按照统一尺寸和可见性配置创建图窗。
- `finiteCdf.m`：忽略非有限值后计算经验CDF坐标，供动态稳定性实验绘图。

### 8.11 `+util`：通用工具

- `makeStream.m`：创建指定基种子和子流编号的独立 `RandStream`，避免算法隐式使用全局随机流。
- `ensureDirectories.m`：创建配置指定的缓存和结果目录。
- `snapshotSignature.m`：计算快拍矩阵的确定性摘要，用于验证所有方法确实使用同一个观测，而不是作为科学指标。

## 9. `tests`：自动化测试

测试按照实现阶段组织，而不是按照单个函数零散组织。

### `TestTask1Contracts.m`

验证配置结构、三种profile、默认参数、任意干扰数、矩形阵列派生维度、输出目录、随机流复现性，以及根入口不允许散落实验数值。

### `TestTask2GeometryScenario.m`

验证旋转矩阵正交性和行列式、零姿态恒等映射、世界—阵列坐标变换、UPA阵元顺序、姿态过程稳态统计、零抖动不变性、快拍维度、QPSK功率以及大快拍SCM向理论协方差收敛。

### `TestTask3BaselinesIncm.m`

验证MVDR和DL-MVDR、指标公式、Oracle上界、二维Capon峰值、可见域掩码、缓存网格流形与逐点参考的一致性、静态INCM的噪声估计、Hermitian性、半正定性和边界场景。

### `TestTask4TemporalRobust.m`

验证二维谱对齐、弱证据回退、精确时序近端解、峰值选择、非重叠自适应求积、结构投影、五项消融开关、离线误差标定和鲁棒约束；同时验证 `epsilonS=0` 且 `gamma>0` 时严格退化为标准MVDR。

### `TestTask5Experiments.m`

验证六方法和六消融变体的配对运行、随机子流独立再生、真值白名单、试次簇bootstrap、500×500统计规模、五个实验入口、完整结果格式、输出路径隔离、拒绝覆盖、事务回滚、备份清理失败和恢复失败等实验系统边界。

### `TestTask6ClassicJitter.m`

验证经典方法选择配置、NSB线性约束、NLMS/RLS训练参考递推、PI参考阵元
约束、ESP子空间投影、WCO最坏情形响应约束、任意方法子集配对运行以及
独立经典方法抖动扫描的结果文件。

### `TestTask7GuLeshem2012Rec.m`

验证Rec式(6)–(7)与手工矩阵计算一致、论文例1参数和角域积分总测度正确、
比较方法只来自论文与现有实现的交集、论文复现同时输出SNR和快拍扫描、UPA
迁移输出网格与扇区诊断，以及非法方法、可能导致原始Rec奇异的快拍数和
不合规二维网格会被配置层拒绝。

## 10. 已实现的比较方法

固定比较顺序为：

1. Oracle MVDR：真实INCM和真实期望导向矢量，仅作上界；
2. 冻结初始MVDR：首周期估计权值，后续保持不变；
3. SCM-MVDR：逐周期直接使用含期望信号的样本协方差；
4. DL-MVDR：逐周期对样本协方差进行配置化对角加载；
5. 静态Capon-INCM：当前周期全局Capon谱、固定期望区排除和均匀积分；
6. 时序稀疏结构约束INCM：谱对齐、时间连续稀疏增强、局部非均匀求积和结构投影。

静态与时序INCM使用同一个鲁棒权值求解接口以及同一组 `gamma`、`epsilonS`，避免把权值求解器差异误认为INCM重构收益。

## 10.1 独立经典方法退化实验

经典方法实验支持以下11种方法标识：

```matlab
["oracle", "nsb", "frozen_mvdr", "scm_mvdr", "nlms", "rls", ...
 "pi", "dl_mvdr", "esp", "wco", "static_incm_mvdr"]
```

可以在调用前自由选择任意非空且不重复的组合，程序保持用户给出的顺序：

```matlab
cfg = config("smoke");
cfg.experiments.classicJitterMethodNames = ...
    ["oracle", "nsb", "scm_mvdr", "dl_mvdr", "wco", ...
     "static_incm_mvdr"];
cfg = uav.config.validate(cfg);
result = exp06_classic_method_jitter_sweep(cfg);
```

NSB使用配置中的标称期望和干扰方向；NLMS与RLS使用已知期望QPSK训练
序列；PI固定配置参考阵元；ESP使用配置的信号子空间维数；WCO采用离线
固定的球形导向不确定集且不附加 `gamma` 正则项；`static_incm_mvdr`
使用静态Capon-INCM重构和标准MVDR，因此不会把鲁棒权值收益混入普通
INCM基线。除Oracle外，所有方法均不能读取当前试次真实姿态、真实方向
或真实INCM。

## 10.2 经典方法参数筛选与锁定验证

`runClassicJitterParameterStudy.m`在同一批配对观测上分别筛选DL-MVDR、
ESP、WCO和静态INCM参数，并在筛选后调用 `runClassicJitterSweep.m`生成
正式实验6结果。所有候选范围集中配置在
`cfg.experiments.classicTuning`，实验函数内部不散落数值参数。

推荐按“smoke筛选、preliminary锁定验证”的顺序运行：

```matlab
startup;
smokeStudy = runClassicJitterParameterStudy(config("smoke"), []);
preliminaryStudy = runClassicJitterParameterStudy( ...
    config("preliminary"), smokeStudy.selections);
```

评分规则固定为“所有抖动扫描点中位SINR的最小值”，优先保证整条曲线
稳定。静态INCM按谱加载、期望区半宽和收缩系数的固定顺序执行坐标筛选；
其他三种方法分别在自身候选范围内取最优值，不能故意选择较差参数。
WCO候选若使离线标定半径达到或超过参考导向矢量范数，将在候选表中明确
标记为不可行，不参与排名，也不会静默删除。

相关模块包括：

- `screenClassicParameters.m`：无绘图候选评价、评分和参数锁定；
- `applyClassicParameterSelection.m`：把smoke锁定参数原样应用到其他profile；
- `runClassicJitterParameterStudy.m`：组织两种快拍数、正式实验和分目录输出。

研究结果分别保存到
`results/classic_parameter_study/<profile>/L8`和 `L64`，总排名表保存在
对应profile目录。当前实测结论和解释边界见
[经典方法参数研究结果](docs/classic_parameter_study_results.md)。

## 10.3 Gu–Leshem 2012 Rec复现与UPA迁移

文献参数集中在 `cfg.literature.guLeshem2012`。`methodNames`必须包含
`gu_leshem2012_rec`，并可在 `oracle`、`esp`、`wco` 中选择对照方法；不得
加入原论文未比较或工程尚未实现的方法。运行方式为：

```matlab
startup;
cfg = config("smoke");
paperResult = exp07_gu_leshem2012_reproduction(cfg);
upaResult = exp08_gu_leshem2012_upa_transfer(cfg);
```

原论文复现与当前UPA迁移的公式边界、参数含义和结果解释规则见
[Gu–Leshem 2012 Rec复现说明](docs/gu_leshem2012_rec_reproduction.md)。

## 11. 结果文件

每个实验成功后生成：

- `results/raw/*.mat`：完整配置、元数据和逐试次逐周期原始结果；
- `results/tables/*.csv`：便于检查和写报告的汇总表；
- `results/figures/*.fig`：MATLAB可编辑图；
- `results/figures/*.pdf`：矢量PDF；
- `results/figures/*.png`：预览图；
- `results/logs/*.json`：结果状态和文件路径清单。

经典参数研究采用独立的 `results/classic_parameter_study` 子目录，避免
覆盖用户直接运行实验6产生的同名结果。

结果发布采用“全部暂存完成后再提交”的方式，避免部分文件成功、部分文件失败却被误认为完整实验。

## 12. 结果解释边界

`smoke`运行通过只能证明模型、算法、测试和结果导出链能够一致执行，不能证明拟研究方法优于基线。

只有当 `preliminary` 或 `paper` 的配对蒙特卡洛结果及置信区间支持差异时，才能写“时序稀疏INCM提高了输出SINR或稳定性”。如果差异不显著，应表述为“完成了模型和初步算法实现，并观察到某种趋势”。所有结论都只适用于当前配置的窄带、单径、普通邻区同频干扰模型及其参数扫描范围。

更严格的公式、维度和实验约定见：

- [模型契约](docs/model_contract.md)
- [变量字典](docs/variable_dictionary.md)
- [实验协议](docs/experiment_protocol.md)
- [Gu–Leshem 2012 Rec复现说明](docs/gu_leshem2012_rec_reproduction.md)
