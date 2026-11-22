# 9×9 UPA原始积分INCM：零抖动理论协方差与Jacobian检查

## 目的与隔离范围

本诊断回答：在当前几何、没有姿态误差和有限快拍误差时，`incm_integral_origin`的二维积分能否通过加密网格收敛？将方向余弦面积测度换成立体角测度能否消除其性能损失？

不研究`incm_integral_processed`，不改变生产算法，也不修改实验6/7/8的已有结果。所有输出写入自动新建的`results/diagnostics/upa_9x9_population_北京时间_随机后缀/`。不同运行不会复用目录。

## 固定场景

- 阵列：9×9 UPA，半波长间距，沿用主工程阵元顺序与幅度归一化。
- 世界方向：期望`[0,-35]`度；两个干扰`[40,-25]`、`[-55,-20]`度。
- 名义姿态与姿态抖动均为零；名义期望导向等于真实期望导向。
- 噪声功率1，两个干扰INR均为30 dB；输入SNR扫描`-10:5:30` dB。
- 保留主工程25度三维球冠期望排除区；所有坐标规则使用同一个连续方向筛选函数。
- 输入为精确理论总协方差，不生成快拍，不做蒙特卡洛，不受`L`影响：

$$
\mathbf R_x=p_s\mathbf a_s\mathbf a_s^H+
\sum_jp_j\mathbf a_j\mathbf a_j^H+\sigma_n^2\mathbf I.
$$

另保存理论INCM和纯噪声协方差两组控制输入，只用于诊断，不冒充可观测算法。

## 保持不变的原始公式

$$
P_C(\mathbf u)=\frac1{\mathbf a^H(\mathbf u)\mathbf R_x^{-1}\mathbf a(\mathbf u)},
\qquad
\widehat{\mathbf R}_{i+n}^{(\mu)}=
\int_{\Omega\setminus\Theta_s}P_C(\mathbf u)\mathbf a(\mathbf u)\mathbf a^H(\mathbf u)\,d\mu.
$$

不对Capon输入加载，不扣谱底，不添加噪声单位阵，不估计功率，不修正期望方向。输出使用同一标准MVDR接口；记录该接口是否触发最小数值加载。

## 坐标与Jacobian

设方位角为$\phi$，高度角为$\alpha\in[-\pi/2,0]$，

$$
\xi_x=\cos\alpha\cos\phi,\quad
\xi_y=\cos\alpha\sin\phi,\quad z=\sin\alpha\le0.
$$

则

$$
d\Omega=\cos\alpha\,d\phi\,d\alpha,
\qquad
d\xi_x\,d\xi_y=|\sin\alpha|\cos\alpha\,d\phi\,d\alpha
=|z|\,d\Omega.
$$

使用$\phi$和$t=-z\in[0,1]$作为求积坐标时：

$$
d\Omega=d\phi\,dt,\qquad
d\xi_x\,d\xi_y=t\,d\phi\,dt.
$$

因此在相同球面节点、相同Capon值上，可以无奇异地比较两种测度，不必在圆盘边缘计算$1/\sqrt{1-\xi_x^2-\xi_y^2}$或人为裁剪它。

必须区分两个问题：同一测度下换坐标必须一致；保持同一个数值$P_C$而更换测度则是另一种重构估计器，不能仅凭Jacobian公式宣称后者就是唯一正确的功率密度模型。

## 网格、收敛标准与独立交叉检查

配置入口：`cfg.experiments.upaPopulationConvergence`。

1. 原生产方向余弦端点网格：每轴`[51,101,201,401,801,1601,3201,6401]`点，使用原始面积权重。
2. 方位—负z中点网格：同样逐级加密；同一批节点分别计算面积积分和立体角积分。
3. 主序列未通过时，按`[9601,12801,16001,19201,25601]`继续追加网格；不得根据结果放宽预设误差阈值。这些追加层级是在6401级的矩阵误差仍未达标后补入，属于增加验证强度，不改变算法与误差标准。
4. 最高精度额外错开半个方位网格，检查峰值与节点偶然对齐的影响。
5. 方位—高度角网格独立计算两级，与方位—负z结果交叉验证。
6. 同一测度最后两次加密都要求：全部9个SNR点最大SINR变化不超过0.05 dB，全部11组输入的最大相对矩阵变化不超过0.001。
7. 坐标交叉检查采用相同0.05 dB与0.001阈值。积分域面积与立体角另对照解析值，避免权重整体错误恰好被MVDR尺度不变性掩盖。

球冠完全位于下半球，本几何下积分补域的解析测度为：

$$
\mu_{xy}=\pi-\pi\sin^2\beta\,|z_s|,
\qquad \mu_{\Omega}=2\pi\cos\beta,\qquad \beta=25^\circ.
$$

分块内核利用矩形UPA的重复阵元差，对同一求和作代数等价加速；不插值谱、不拟合功率。正式运行前以9×9高SNR输入逐元素核对直接原式；单元测试还检查白噪声下的球面协方差解析表达。

## 运行及输出

```matlab
startup;
cfg = config('paper');
result = diagnose_upa_population_convergence(cfg);
result.outputDirectory
result.allChecksPassed
```

输出包括`inputs.mat`、`checkpoint.mat`、`study.mat`、`parity.csv`、`summary.csv`、`convergence_checks.csv`、`coordinate_checks.csv`、`final_comparison.csv`、北京时间`progress.log`及两组FIG/PDF/PNG图。`allChecksPassed=false`表示尚未完成预设标准下的数值收敛，不能只根据最高精度单条曲线宣布收敛。

需要续算时，把`cfg.experiments.upaPopulationConvergence.resumeCheckpoint`设为已完成运行的`checkpoint.mat`绝对路径。程序核对理论输入与排除区完全一致，再新建目录接续；不会修改原检查点。临时性能设置不改变求积节点、公式或误差阈值。

### 可选本地加速

默认`nativeAcceleration=false`，不要求MEX或额外编译器。此次大网格续算采用`upaPopulationQuadratureNative.cpp`的OpenMP并行双精度求和。9×9、全部11组输入、51及401级网格下，四种规则与MATLAB的相对矩阵差不超过`2.2e-12`，纳入节点数一致。生产重构算法不调用该内核。

当前机器的可复现编译命令（PowerShell，在工程根目录执行；输出目录必须先新建）为：

```powershell
& 'C:/msys64/mingw64/bin/g++.exe' -O3 -march=native -ffp-contract=off -fopenmp -shared -DMATLAB_MEX_FILE -I D:/MATLAB/extern/include src/+uav/+incm/upaPopulationQuadratureNative.cpp -L D:/MATLAB/extern/lib/win64/microsoft -llibmex -llibmx -o results/diagnostics/native_build_20261114/upaPopulationQuadratureNative.mexw64
```

`-ffp-contract=off`很重要：避免可见圆盘边界的乘加融合改变个别边界点的纳入判定。不能使用会改变浮点语义的`-ffast-math`。

```matlab
addpath(fullfile(pwd,'results','diagnostics','native_build_20261114'));
cfg.experiments.upaPopulationConvergence.nativeAcceleration = true;
cfg.experiments.upaPopulationConvergence.nativeThreadCount = 16;
```

记录中的每个网格均保存是否使用本地加速；续算前已经完成的网格仍保留原记录及运行时间。

## 结果与结论

最终交付目录为：

`D:/workspace/UAVjitter/results/diagnostics/upa_9x9_population_20261114_101512_283_tpce39776d_9a12_4fcd_a077_cea300ff535a/`

该目录复用同场景检查点完成最后一次校验、真值控制表及绘图导出，不重跑已完成求积。检查点来源链记录在配置和元数据中。源目录的已完成网格及数据均保留，未覆盖；原有实验8的两份MAT与两份UPA/ULA FIG经SHA-256核对完全不变。

这里的确定性曲线不使用中位数、置信区间或误差棒；其横纵值均来自同一固定场景的理论协方差。5项针对性单元/流程测试通过，最终15项数值检查全部通过，`result.allChecksPassed=true`。未运行与本诊断无关的全量蒙特卡洛实验。

### 1. 网格已经收敛，差距并未消失

SNR=10 dB时，原始方向余弦面积积分的输出如下；这里每轴点数为G，总网格为G×G，再筛选可见域及排除区。

| G | 输出SINR（dB） |
|---:|---:|
| 51 | 22.333925 |
| 101 | 21.679118 |
| 201 | 21.949919 |
| 401 | 21.912739 |
| 801 | 21.908558 |
| 1601 | 21.932808 |
| 3201 | 21.927239 |
| 6401 | 21.929769 |
| 9601 | 21.930160 |
| 12801 | 21.929857 |
| 19201 | 21.930324 |

原始网格从12801到19201，全部SNR点的最大SINR变化仅`0.0004663 dB`；最大相对矩阵变化`6.576e-6`。从201到最终网格，SNR=10 dB的SINR只改变约`0.0196 dB`。因此，在本场景零抖动、理论协方差条件下，不能再把约7 dB的性能损失主要归因于201网格不够密。

注意51网格反而更高不表示其更准确；粗网格结果非单调，可能漏算或偶然踩中尖峰，必须与最终收敛值比较。

### 2. Jacobian不能消除剩余损失

三组积分最终均使用19201级网格。以SNR=10 dB为例：

| 重构/权值输入 | 输出SINR（dB） | 与Oracle差距（dB） |
|---|---:|---:|
| 真实INCM＋准确期望导向（Oracle） | 29.059770 | 0 |
| 原始方向余弦面积积分 | 21.930324 | 7.129446 |
| 相同球面节点上的面积积分 | 21.930399 | 7.129371 |
| 相同球面节点上的立体角积分 | 22.781572 | 6.278197 |

立体角积分相对同节点面积积分改善`0.851173 dB`。全体SNR点的改善范围为`0.849283–0.851195 dB`，并没有接近Oracle。

如果定义目标为立体角积分，换到方位—高度角时必须乘`cos(高度角)`，换到方向余弦时必须带`1/|z|`。但当前生产实现定义的本来就是方向余弦面积积分，它不需要凭空再乘一个立体角Jacobian；保持同样的Capon数值而更换测度，是改变了重构估计器，不是单纯把一个错误的坐标换元补对。

本次实验证明：**正确的立体角权重可以带来有限改善，但“遗漏Jacobian”不是本场景性能差距的主要根因。**

### 3. 交叉检查排除了什么

- 同一面积测度，方向余弦坐标与方位—负z坐标：全部SNR点最大差`7.52e-5 dB`，最大相对矩阵差`1.67e-5`。
- 同一立体角测度，方位—负z与方位—高度角坐标：最大差`7.39e-6 dB`，最大相对矩阵差`7.98e-7`。
- 19201级网格错开半个方位单元：最大SINR差不超过`2.58e-5 dB`。
- 独立方位—高度角9601到19201加密：最差相对矩阵变化`0.0007843`，小于预先设置的`0.001`。
- 面积与立体角补域的解析测度相对误差分别约`9.91e-8`、`4.14e-8`。
- 所有输出权值均未触发MVDR最小数值加载。

因此，当前结果不是网格偶然正好命中某个干扰峰，也不是方位/高度角换元权重写错。解析域检查的SINR字段为`NaN`表示该检查只比较面积、SINR不适用；不是算法输出了非有限SINR。

### 4. 两个只用于排查的真值对照

详见`truth_controls.csv`；不可将这些使用真值的控制组冒充实际算法。

1. 完全移除期望信号，给同一积分公式输入真实INCM：SNR=10 dB时，面积积分`21.971947 dB`、立体角积分`22.806323 dB`。相比总协方差输入仅增加约`0.042/0.025 dB`，损失仍明显存在。所以本轮准确导向条件下，“输入总协方差含期望信号”不是主要解释。
2. 使用真实干扰方向读取Capon点功率，按离散外积求和并加入已知真实噪声单位阵：SINR与Oracle在数值精度内一致，点功率相对真实干扰功率的误差约`1.24e-5`。此控制同时使用真实方向和噪声功率，不能解释为实际峰值搜索方法已经达到Oracle；它只说明UPA导向、Capon点值及SINR接口本身能给出正确结果。

综合判断：在当前9×9、固定几何和25度排除区下，问题更指向“把有限阵列Capon谱直接当连续功率密度作外积积分”所得协方差与真实离散源INCM的差异，而不是简单的网格不够密或漏乘Jacobian。本次尚未分离峰面积标定、谱底/旁瓣与排除区的各自贡献，不能把其中某一项宣布为唯一根因，也不能据此声称所有UPA角域积分方法都不可行。

### 5. 产物定位与图形核对

- `grid_convergence.fig/.pdf/.png`：SNR=10 dB的网格收敛；横轴为对数坐标，不同坐标和测度用不同线形/标记。
- `sinr_comparison.fig/.pdf/.png`：全体SNR点最高网格结果。
- `summary.csv`：每网格、每测度、每SNR的SINR、泄漏、噪声输出、加载和矩阵形状距离。
- `convergence_checks.csv`、`coordinate_checks.csv`：全部15项数值检查。
- `truth_controls.csv`：真值辅助的原理控制，不参与实际方法排名。
- `study.mat`：包含完整配置、所有网格重构矩阵、理论输入、指标和元数据。
- `targeted_tests.mat`：5项针对性测试结果。

已逐图检查中文标签、图例、对数横轴、曲线排序和FIG/PDF/PNG导出。面积积分的不同坐标曲线几乎重合是换元一致性的数值证据，不是遗漏曲线。PDF为矢量输出。诊断是固定理论输入，无随机重复，因而不画误差棒。
