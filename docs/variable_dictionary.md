# 变量字典

| 符号或字段 | 尺寸 | 含义 | 可访问范围 |
| --- | ---: | --- | --- |
| `Nx, Ny` | 标量 | UPA沿阵列坐标系x、y方向的阵元数 | 配置参数 |
| `M` | 标量 | 阵元总数 `Nx*Ny` | 派生配置 |
| `J` | 标量 | 根据干扰方向和INR数组推导的干扰数 | 派生配置 |
| `K, L` | 标量 | 处理周期数和每周期快拍数 | 配置参数或扫描参数 |
| `eta, DeltaEta` | `3×K` | `[横滚角;俯仰角;偏航角]`姿态及其偏差，单位为弧度 | 场景真值 |
| `R_k` | `3×3` | 从阵列坐标系到世界坐标系的 `Rz Ry Rx` 旋转矩阵 | 场景真值 |
| `uWorld` | `3×(J+1)` | 世界坐标系中的固定单位方向 | 配置先验 |
| `uArray_k` | `3×(J+1)` | 当前周期在阵列坐标系中表示的方向 | 真值或诊断量 |
| `xi_k` | `2×(J+1)` | 阵列坐标系x、y方向余弦 | 真值或诊断量 |
| `A_k` | `M×(J+1)` | 真实导向矢量矩阵，第一列为期望信号 | Oracle或指标模块 |
| `aRef` | `M×1` | 名义期望导向矢量参考 | 算法先验 |
| `X_k` | `M×L` | 复数阵列接收快拍 | 算法观测 |
| `Rx_k` | `M×M` | 样本协方差矩阵 `X X^H/L` | 算法观测 |
| `Rin_k` | `M×M` | 理论干扰加噪声协方差矩阵 | Oracle或指标模块 |
| `RinHat_k` | `M×M` | 重构得到的INCM | 算法估计量 |
| `incm.staticCovfit.loadingRatio` | 标量 | `static_incm_mvdr`候选Capon谱的相对加载系数，当前默认0.12 | 配置参数 |
| `incm.staticCovfit.fitRegularization` | 标量 | 非负协方差拟合的L2正则化系数，当前默认5000 | 配置参数 |
| `incm.staticCovfit.peakThresholdRatio` | 标量 | 相对谱峰阈值，当前默认0.15 | 配置参数 |
| `incm.staticCovfit.peakNmsRadiusBins` | 非负整数 | 二维非极大值抑制半径，单位为网格点，当前默认0 | 配置参数 |
| `incm.staticCovfit.maximumCandidateCount` | 正整数 | 干扰候选数配置上限，实际还受 `M-1` 限制，当前默认7 | 配置参数 |
| `p_k, c_k` | 网格矩阵 | 残余Capon谱和时序增强谱 | 估计量或状态量 |
| `G_k` | 平移诊断量 | 从上一周期空间谱到当前周期空间谱的对齐关系 | 估计量或状态量 |
| `w` | `M×1` | 复波束形成权值 | 算法输出 |
| `sinrLinear` | 标量或数组 | 期望输出功率与干扰加噪声输出功率之比 | 主要原始指标 |
| `sinrDb` | 标量或数组 | 在线性域汇总后执行 `10*log10` 得到的指标 | 报告指标 |
| `desiredGain` | 标量或数组 | 期望方向增益 `|w^H a_s|^2` | 指标 |
| `interferenceLeakage` | 标量或数组 | 不含接收机噪声的总干扰输出功率 | 指标 |
| `desiredSteeringMismatchNormalized` | 标量或数组 | 式（3-10）定义的期望信号归一化导向矢量失配 | 场景评价指标 |
| `incmNmse` | 标量或数组 | 相对Frobenius误差的平方 | 指标 |
| `exp02MethodNames` | `1×N_method`字符串 | 实验2实际计算、绘图和导出的方法子集及顺序 | 配置参数 |
| `classicJitterMethodNames` | `1×N_method`字符串 | 独立经典方法抖动实验实际执行的方法及绘图顺序 | 配置参数 |
| `classicTuning.snapshotCounts` | 正整数向量 | 经典参数研究分别评价的每周期快拍数，当前为8和64 | 配置参数 |
| `classicTuning.scoreRule` | 字符串 | 参数筛选评分规则，固定为各抖动点中位SINR的最小值 | 配置参数 |
| `minimumPointMedianSinrDb` | 标量 | 一条抖动扫描曲线中最低的逐点中位SINR | 参数筛选指标 |
| `validCandidate` | 逻辑标量 | 参数候选是否满足标定和求解可行性 | 参数筛选诊断 |
| `selection` | 结构体 | 替换前参数研究对DL、ESP、WCO和旧静态Capon积分INCM锁定的参数；其中旧静态INCM字段不适用于当前CovFit | 历史smoke输出/验证输入 |
| `nlmsState` | 结构体 | NLMS上一周期权值 | 算法状态 |
| `rlsState` | 结构体 | RLS上一周期权值和逆相关矩阵 | 算法状态 |
| `epsilonS` | 标量 | WCO或鲁棒权值使用的球形导向不确定集半径 | 离线标定配置 |

信号源列顺序固定为期望信号在前，随后按配置顺序排列各干扰信号。传感器顺序为 `vx kron vy`。实验原始数组的前两维固定为 `试次×周期`，若存在扫描点和方法或消融变体，则依次放在后续维度。同一观测对应的快拍签名会复制到所有比较方法或变体，用于证明配对实验确实使用相同数据。只有不适用于某种方法的INCM指标允许使用NaN；主要SINR、期望增益和干扰泄漏必须为有限值。
