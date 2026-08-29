# 姿态抖动下经典方法退化实验实施计划

> **执行要求：** 使用 `superpowers:executing-plans` 在当前会话逐项实施；每项均遵循测试先行、最小实现和回归验证。

**目标：** 新增一个可自由选择 Oracle、NSB、两种传统MVDR、NLMS、RLS、PI、DL、ESP、WCO和普通静态INCM-MVDR的独立姿态抖动扫描实验。

**架构：** 保留现有五组实验及固定六方法合同不变。新增经典方法专用配置、波束形成函数和实验运行器；所选方法列表仅影响新实验，并按配置顺序动态生成数据维度、图例和结果表。

**技术栈：** MATLAB R2024a、`matlab.unittest`、Optimization Toolbox、现有 `+uav` 包结构。

**依据：** 开题报告表2-1及第2.1.1节；当前用户批准的方法集合与“原始经典方法、不引入改进”的约束。

## 全局约束

- 不修改 `D:/workspace/hithesis-open` 中的开题报告源文件。
- 不覆盖用户现有 `config.m` 和 `runExp02.m` 修改。
- 非Oracle方法只能读取 `algorithmObservation` 白名单、配置先验和自身历史状态。
- 所有数值参数进入 `config.m`，算法文件不得散落实验常量。
- 所有新增文档、函数头、测试意图和行内说明使用中文。
- 新实验不加入 `run_all`，避免默认批处理突然增加高开销实验。

---

### 任务1：经典方法选择配置与合同

**文件：**
- 修改：`config.m`
- 修改：`src/+uav/+config/validate.m`
- 创建：`tests/TestTask6ClassicJitter.m`

**接口：**
- 产生：`cfg.experiments.classicJitterMethodNames`，允许从固定合法集合中选择任意非空且不重复的组合。
- 独立实验复用：`cfg.experiments.periodCounts.jitterSweep`，避免改变原有周期配置契约。
- 产生：`cfg.beamforming.classic.nlms`、`rls`、`pi`、`esp` 参数。

- [x] 写配置字段、任意合法子集、重复名称和未知名称的失败测试。
- [x] 运行测试，确认因字段缺失而失败。
- [x] 在 `config.m` 增加默认11方法列表及 NLMS 步长/正则项、RLS 遗忘因子/初始化对角量、PI参考阵元、ESP子空间维数。
- [x] 在 `validate.m` 校验非空唯一子集、周期数和经典参数范围。
- [x] 运行配置测试并确认通过。

### 任务2：经典波束形成模块

**文件：**
- 创建：`src/+uav/+beamforming/nullSteering.m`
- 创建：`src/+uav/+beamforming/adaptiveReference.m`
- 创建：`src/+uav/+beamforming/powerInversion.m`
- 创建：`src/+uav/+beamforming/eigenspaceProjection.m`
- 修改：`src/+uav/+beamforming/solve.m`
- 修改：`tests/TestTask6ClassicJitter.m`

**接口：**
- `nullSteering(nominalSteering,cfg)`：返回满足标称期望单位响应及标称干扰零响应的最小范数NSB权值。
- `adaptiveReference(method,obs,priorState,cfg)`：按快拍执行经典NLMS或RLS，参考信号为 `sqrt(desiredPower)*pilot`。
- `powerInversion(R,cfg)`：固定配置参考阵元响应为1并最小化输出功率。
- `eigenspaceProjection(R,aRef,cfg)`：将标称导向矢量投影至SCM主特征子空间后求标准MVDR。
- `solve("wco",R,aRef,cfg)`：使用球形不确定集和零附加正则项求经典最坏情形权值。

- [x] 写NSB标称约束、NLMS/RLS确定性收敛、PI参考阵元约束、ESP投影空间及WCO响应约束测试。
- [x] 运行测试，确认函数不存在或方法不支持导致预期失败。
- [x] 分别实现五个经典模块，数值修复只用于避免奇异矩阵产生NaN，不改变经典目标。
- [x] 运行模块测试并确认全部通过。

### 任务3：可选组合运行器与独立实验

**文件：**
- 创建：`src/+uav/+experiment/runClassicMethods.m`
- 创建：`src/+uav/+experiment/collectClassicMethods.m`
- 创建：`src/+uav/+experiment/runClassicJitterSweep.m`
- 创建：`experiments/exp06_classic_method_jitter_sweep.m`
- 修改：`tests/TestTask6ClassicJitter.m`

**接口：**
- `[record,nextState]=runClassicMethods(obs,priorState,cfg,methodNames)`。
- `[raw,seeds]=collectClassicMethods(cfg,experimentIndex,pointIndex,periodCount,methodNames)`。
- `result=exp06_classic_method_jitter_sweep(cfg)`。

- [x] 写任意方法子集保持配置顺序、共享快拍摘要、非Oracle真值隔离和普通INCM使用标准MVDR的失败测试。
- [x] 运行测试并确认新入口缺失导致失败。
- [x] 实现方法调度：Oracle、NSB、冻结MVDR、SCM-MVDR、NLMS、RLS、PI、DL-MVDR、ESP、WCO、静态Capon-INCM-MVDR。
- [x] 实现配对蒙特卡洛扫描、线性域统计、SINR中位数/第5百分位双面板、MAT/CSV/FIG/PDF/PNG/JSON导出。
- [x] 运行新测试类并确认通过。

### 任务4：中文说明与完整验证

**文件：**
- 修改：`README.md`
- 修改：`docs/variable_dictionary.md`
- 修改：`docs/experiment_protocol.md`

- [x] 说明11种经典方法、配置选择示例、各方法先验和结果解释边界。
- [x] 对新增和修改MATLAB文件运行 `checkcode`。
- [x] 运行 `runtests('tests/TestTask6ClassicJitter.m')`。
- [x] 运行完整测试套件并区分本次失败与既有配置/测试不一致。
- [x] 用紧凑配置运行一次独立实验，确认所选子集与六类结果文件完整生成。
