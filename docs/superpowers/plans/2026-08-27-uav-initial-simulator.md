# UAV姿态抖动下鲁棒INCM仿真系统实施计划

> **面向自动化开发执行者：** 每项任务必须严格遵循“红—绿—重构”的测试驱动开发流程。不得分派子代理。每完成一项任务后单独提交。

**目标：** 基于MATLAB R2024a建立可复现的仿真系统，用于研究UAV姿态抖动引起的阵列失配、同频干扰、基线波束形成、静态Capon-INCM重构，以及时序稀疏和结构约束INCM重构。

**架构：** `config.m` 是唯一数值参数来源。`src/+uav` 下的纯包函数分别实现几何、场景、算法、指标和绘图；实验入口只负责组织配置中规定的参数扫描。真值与算法输入严格隔离，只允许Oracle方法和指标模块读取。

**技术环境：** MATLAB R2024a、MATLAB Test、Optimization Toolbox，可选Parallel Computing Toolbox。不使用CVX，不生成5G完整波形，不包含学习子系统。

**方案依据：** 用户在对话中批准的方案，并以 `D:/workspace/hithesis-open/examples/hitart/reports/body/report_shenzhen_doctor_opening.tex` 和 `D:/workspace/hithesis-open/第三章研究内容及方案计划.md` 为依据。

## 全局约束

- 只能在 `D:/workspace/UAVjitter` 中工作，不得修改HiThesis开题报告工程。
- 默认主阵列为 `3×3`，但所有正式函数必须支持 `Nx~=Ny` 的矩形阵列。
- 干扰数 `J` 根据配置的方向和功率数组推导；不得写死循环或矩阵尺寸，并应支持 `J=0` 和任意 `J`。
- 物理模型采用窄带、远场、单径和世界坐标系信号源固定的假设；只有UAV姿态变化会改变信号在阵列坐标系中的表观方向。
- 当前报告场景是正常邻区基站同频干扰，不是压制式干扰。
- 所有运行数值和实验参数都必须放在 `config.m` 中；测试可以使用显式夹具参数。
- 角度只能在配置边界使用度，模型内部统一使用弧度。
- 真值不得进入非Oracle算法函数。
- MATLAB代码必须显式传递 `RandStream` 以保证可复现性，并避免使用全局随机数状态。
- 不包含学习模型、5G NR完整波形、多径、信号源运动、传感器输入，也不修改开题报告源文件。

---

### 任务1：工程骨架、配置与契约

**涉及文件：** 创建 `config.m`、`startup.m`、`run_all.m`、`README.md`、`src/+uav/+config/validate.m`、`src/+uav/+util/makeStream.m`、契约文档、输出和缓存目录，以及针对性的测试。

**要求：**

- `config(profile)` 只接受 `smoke`、`preliminary` 或 `paper`，并返回以下嵌套配置分组：`random`、`array`、`timeline`、`attitude`、`sources`、`signal`、`spectrum`、`incm`、`beamforming`、`monteCarlo`、`experiments`、`plot` 和 `output`。
- 默认参数：`Nx=Ny=3`，半波长阵元间距；期望方位角/高度角 `[0,-35]` 度；干扰方向 `[40,-25;-55,-20]` 度及INR `[20,15]` dB；SNR为10 dB；处理周期1 ms；`L=8`；相关系数 `[0.95,0.95,0.98]`；基种子20260827；抖动扫描 `[0,0.5,1,2,4,8]` 度；快拍扫描 `[2,4,8,16,32,64]`。
- 三种运行规模只能改变工作量和网格分辨率：`smoke` 为3次试验及较短的 `K`，`preliminary` 为100次试验，`paper` 为500次试验。
- 配置校验需要推导 `M`、`J`、信号功率、世界坐标系单位方向、名义导向矢量参考和输出路径；同时校验向量长度、有限数值、可见下半球信号源方向和正尺寸。
- `startup` 只添加 `src`、`experiments` 和 `tests` 路径。
- `run_all(profile)` 先运行测试，再依次运行实验。
- 在API尚未实现时测试必须先失败，实现后再通过。

### 任务2：几何、姿态、矩形UPA与快拍生成

**涉及文件：** 在 `+geometry`、`+array` 和 `+scenario` 下创建职责集中的包函数，并增加测试。

**接口：**

- `R = uav.geometry.rotationZYX(etaRad)`。
- `uArray = uav.geometry.worldToArray(R,uWorld)`；`[xi,azElRad] = uav.geometry.directionCoordinates(uArray)`。
- `a = uav.array.steeringUPA(xi,Nx,Ny,dOverLambda)`，传感器顺序为 `vx kron vy`，阵元幅度为1，导向矢量范数为 `sqrt(M)`。
- `scenario = uav.scenario.generate(cfg,stream)` 返回稳态AR(1)姿态偏差、旋转矩阵、方向和真实导向矢量，不预生成全部快拍。
- `obs = uav.scenario.generatePeriod(cfg,scenario,k,stream)` 返回 `X`、`Rx`、已知QPSK训练序列、算法安全字段，以及独立的真值结构；真值结构包含真实导向矢量、真实INCM、信号功率和方向。

**要求：** 使用稳态关系 `Qnu=Qeta-F*Qeta*F'`；期望信号和干扰QPSK序列相互独立；噪声为复白噪声；不引入衰落；支持矩形阵列以及 `J=0`、`J>=M`。测试应覆盖旋转矩阵正交性、零姿态、导向矢量字面值、稳态协方差、大快拍协方差收敛、随机流确定性、`L<M` 以及无NaN/Inf。

### 任务3：基线、指标、Capon空间谱与静态INCM

**涉及文件：** 在 `+beamforming`、`+metrics`、`+spectrum` 和 `+incm` 下增加职责集中的函数及测试。

**接口：**

- `[w,diag] = uav.beamforming.solve(method,R,aRef,cfg)`，本任务先支持 `mvdr` 和 `dlmvdr`。
- `metrics = uav.metrics.evaluate(w,truth,cfg)` 返回线性域和dB域输出SINR、期望增益、干扰泄漏，以及可选的INCM归一化均方误差。
- `[powerMap,diag] = uav.spectrum.loadedCapon(Rx,cfg)` 在配置的下半球方向余弦网格上扫描。
- `[Rhat,state,diag] = uav.incm.reconstruct(method,obs,priorState,cfg)`，本任务先支持 `static`。

**要求：** 实现Oracle MVDR、冻结初始MVDR的调度支持、SCM-MVDR、DL-MVDR、固定名义期望信号排除区域、均匀求积、较小特征值噪声估计，以及Hermitian、半正定和收缩结构约束。测试覆盖Capon谱峰定位、半正定和Hermitian输出、Oracle SINR上界、`J=0` 以及奇异SCM加载。

### 任务4：时序稀疏INCM、空间谱对齐、自适应求积与鲁棒权值

**涉及文件：** 增加时序空间谱、谱对齐、谱峰及局部细化、结构投影和鲁棒权值函数，并增加测试。

**要求：**

- 使用归一化相关估计有界整数二维空间谱平移；证据较弱时使用恒等映射。
- 在 `c>=0` 约束下，采用精确的逐元素分段近端解求解 `||p-c||^2 + lambdaS||c||1 + lambdaT||c-GcPrev||1`。
- 检测谱峰时不得使用真实 `J`，谱峰数量上限为 `M-1`；构造互不重叠的“区域外粗节点、区域内细节点”求积节点及面积权重。
- 重构INCM后依次执行Hermitian化、特征值下限修复和收缩。
- 训练序列正交投影只能作为消融开关，默认关闭。
- INCM方法共用的鲁棒权值采用实变量形式，通过Optimization Toolbox求解 `Re(w^H aRef)-epsilonS||w||>=1`；`epsilonS=0` 时必须返回标准MVDR。
- 在 `uav.incm.reconstruct` 中增加 `temporal` 支持，并在波束形成求解器中增加 `robust` 支持。
- 测试应覆盖已知合成平移、弱证据恒等映射、近端解与手工标量算例及数值搜索的比较、无重叠求积、半正定输出、鲁棒约束可行性以及 `epsilonS=0` 回归。

### 任务5：实验、结果导出、文档与端到端验证

**涉及文件：** 创建五个实验函数、共享调度、汇总和绘图辅助函数、文档以及集成测试。

**要求：**

- 实验1：姿态、方向和导向矢量失配时序，以及零抖动验证。
- 实验2：配对抖动扫描，比较Oracle、冻结、SCM、DL、静态INCM和时序INCM；输出SINR中位数、第5百分位、95%置信区间、期望增益和干扰泄漏。
- 实验3：快拍扫描，输出SINR、INCM归一化均方误差、最小特征值、条件数和运行时间。
- 实验4：输出动态SINR时序和累积分布函数、低分位、低于阈值的概率、干扰泄漏、期望增益和运行时间。
- 实验5：依次执行 `lambdaT=0`、恒等谱对齐、均匀求积、关闭结构增强和可选训练序列投影的配对消融。
- MAT文件保存完整配置、MATLAB与工具箱元数据、方法名称、随机种子和子流、逐周期逐试次原始指标、时间戳及运行时间；同时导出CSV汇总表、FIG、矢量PDF和PNG图。
- `smoke` 运行规模必须完成全部五个实验并生成输出，生成结果不得纳入Git版本控制。
- README必须写明准确的测试及 `smoke` 命令，并说明论文结论的证据边界。

**最终验证命令：**

`matlab -batch "startup; r=runtests('tests','IncludeSubfolders',true); assertSuccess(r)"`

`matlab -batch "startup; run_all('smoke')"`
