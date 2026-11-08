# UPA积分迁移的变量、测度、积分域与性能根因

## 一、本轮结论与边界

研究仅在 `codex/upa-integral-jacobian-study` 进行，未修改或合并main；不再研究 `incm_integral_processed`。
正式 `incm_integral_origin` 入口保持原定义，本轮增加独立诊断而非悄悄替换算法。

需要纠正的表述是：**不能把“UPA效果差”直接归结为二维积分不可行，也不能把补Jacobian称为原论文唯一正确的UPA迁移。**
当前证据支持：数值求积已收敛；方向换元实现一致；但把有限孔径Capon谱直接当作连续角功率密度，再用广角域积分重构离散源与独立接收机噪声，存在估计偏差。排除区选择也显著影响本例结果。

本轮主控制组为当前3×3、半波长UPA和原名义几何：期望[0,-35]度，干扰[40,-25]与[-55,-20]度，INR各30 dB，噪声方差1。
使用零姿态抖动的**总体协方差**排除有限快拍和方向失配；这是根因控制，不是把正式实验8的5度RMS改为0，也不是新的蒙特卡洛性能曲线。
全部方法仍用同一名义期望导向与标准MVDR接口。

## 二、到底是不是同时对方位角和高度角积分

是二维方向域积分，不是对三维位置或距离作体积分。定义方位角为 $\phi$、高度角为 $\alpha$：

$$
\mathbf u(\phi,\alpha)=
[\cos\alpha\cos\phi,\ \cos\alpha\sin\phi,\ \sin\alpha]^T.
$$

本模型阵面为x-y平面，方向余弦为 $u=\cos\alpha\cos\phi$、$v=\cos\alpha\sin\phi$，负半球 $z=\sin\alpha\le0$。

$$
d\Omega=\cos\alpha\,d\phi\,d\alpha,
\qquad du\,dv=|\sin\alpha|\cos\alpha\,d\phi\,d\alpha=|z|d\Omega.
$$

因此 $d\Omega=du\,dv/\sqrt{1-u^2-v^2}$，也可用无边缘奇异的坐标写成 $d\Omega=d\phi\,dz$。
若改用从正z轴量起的天顶角 $\theta$，才写作 $d\Omega=\sin\theta\,d\theta\,d\phi$；不能将高度角和天顶角混用。

但权重不能脱离“密度是对什么变量定义的”单独讨论。相同物理功率测度满足

$$
d\mu=S_\Omega\,d\Omega
=S_{\phi\alpha}\,d\phi\,d\alpha
=S_{uv}\,du\,dv,
\quad S_{\phi\alpha}=S_\Omega\cos\alpha,
\quad S_{uv}=S_\Omega/|z|.
$$

只换坐标并同时转换密度，结果应相同；**保持同一个数值Capon谱，仅把权重由du dv换成dOmega或dphi dalpha，是更换估计器，不是同一物理模型的纯换元。**
平坦角度乘积积分不一定“数学错误”，但需要说明相应密度/角度先验，不能因为某条曲线较高就认定其物理定义唯一正确。

## 三、积分域与排除区的合理性审查

当前负半球取 $\phi\in[-\pi,\pi)$、$\alpha\in[-\pi/2,0]$，等价于方向余弦单位圆盘。
这是已知源在阵面下侧的先验，不是所有UPA天然只有半球。
理想无方向性平面阵满足 $\mathbf a(u,v,z)=\mathbf a(u,v,-z)$，仅靠该阵列无法区分镜像来向。
若对整个球面积分并将排除区同步镜像，因响应镜像相同，结果只整体乘2，MVDR权值不变；只排除下半球球冠却保留上半球镜像期望峰，则是不一致的排除。
姿态范围足够大导致源跨越阵面时，需要重新审查“源侧已知”以及镜像排除；本轮零抖动控制不存在跨面问题。

当前排除区是三维球冠而非方位/高度矩形：

$$
\Theta_s=\{\mathbf u:\arccos(\mathbf u^T\bar{\mathbf u}_s)\le\beta\},
\qquad \mathcal D=\{z\le0\}\setminus\Theta_s.
$$

球冠是各向同性方向不确定性的合理几何选择，但**25度并非原论文自动给出的UPA参数**。
Gu–Leshem 2012正文式(7)是一维角域积分，排除区宽度须结合阵列分辨率和环境选择；没有给出当前二维球冠的唯一迁移公式。
当前干扰距名义期望方向分别约35.795度、50.381度，虽然真干扰点均在25度球冠外，也不保证其有限孔径Capon谱周围的全部有效质量不受截除。

在SNR10 dB、总体协方差、零抖动下，立体角积分的排除半径控制结果为：

| 球冠半径 | 输入总协方差的SINR/dB | 输入真实INCM的SINR/dB |
|---|---:|---:|
| 0度 | 15.1903 | 12.4209 |
| 5度 | 15.2196 | 12.4182 |
| 15度 | 14.6823 | 12.1972 |
| 25度 | 9.4597 | 10.4920 |
| 35度 | 0.6951 | 3.7500 |

这里0度不是有效的鲁棒期望排除设计，而是“是否排除”控制；不能据此把正式算法改为不排除期望信号。
无方向失配时，完整精确总协方差也可以生成最优MVDR，因此总协方差积分偶尔优于真实INCM积分不矛盾；积分映射本身是非线性的，不保留这种理论等价性。
5度只是诊断点，不是本轮选定的新默认值，更不能声称它覆盖正式三轴5度RMS的全部方向误差。

## 四、网格到底有多密，还能不能靠加密解决

上一轮最高6401×6401个 $(\phi,z)$ 中点节点，排除前约4097万个节点。
方位步长约0.056241度，z步长约0.000156226；高度角非均匀，不能误称为0.056241度的均匀高度步长。

本轮增加两种独立验证：

- 均匀方位/高度角网格：每维401、801、1601、3201、6401节点，同时计算三种测度。
- 方位/z网格进一步升到12801×12801：方位步长约0.028123度，z步长约0.000078119。

五个SNR总体场景的检查：

| 检查 | 最大SINR绝对变化/dB | 最大矩阵相对差 |
|---|---:|---:|
| 6401方位/高度角 vs 6401方位/z，立体角测度 | 0.0001181 | 0.000001301 |
| 6401 vs 12801方位/z，立体角测度 | 0.0002655 | 0.000002567 |

所以在已检查场景中，约10 dB差距不是尚未分辨窄峰导致的残留求积误差。12801只计算本轮少量总体控制输入，不能宣称上一轮905个输入全部都算到了12801。

以SNR10 dB为例：

| 估计器或控制 | 输出SINR/dB |
|---|---:|
| 原du dv积分，12801 | 7.9477 |
| dOmega积分，12801 | 9.4601 |
| dphi dalpha平坦角度乘积积分，6401 | 10.4831 |
| 真实INCM先输入Capon再作dOmega积分，25度排除 | 10.4922 |
| 在真实干扰方向直接读Capon峰功率再作离散源和，Oracle控制 | 19.3608 |
| 真实INCM直接求MVDR | 19.3608 |

真实INCM输入控制仍失效，说明“SCM有期望信号”“有限快拍”“姿态抖动”都不是这一控制组损失的必要条件。
同时，准确方向处的峰值功率已足以近乎恢复上界；把峰变成连续面积后才出现显著偏差，这是后续应研究的估计映射问题。

## 五、关键原理：Capon峰高不等于角功率密度的质量

物理协方差的安全起点是非负功率测度，而不是先规定必须使用某个坐标网格：

$$
\mathbf R_{i+n}=\int\mathbf a(\Omega)\mathbf a^H(\Omega)\,d\mu_i(\Omega)+\sigma_n^2\mathbf I.
$$

单径离散源对应 $d\mu_i=\sum_j p_j\delta_{\Omega_j}$，自然得到源协方差外积求和。若有连续散射角密度，才转为相应测度下的密度积分。
两者都允许UPA；问题不是二维积分不能表示协方差，而是实际拿什么来估计 $d\mu_i$。

以下是本工程的独立推导，不冒称为原论文的定理。单个源加接收机白噪声时

$$
\mathbf R=\sigma_n^2\mathbf I+p\mathbf a_0\mathbf a_0^H,
\quad \|\mathbf a\|^2=M,
\quad P_C(\mathbf a)=\frac{\sigma_n^2}{M-\frac{p}{\sigma_n^2+pM}|\mathbf a^H\mathbf a_0|^2}.
$$

准确方向处有 $P_C(\mathbf a_0)=p+\sigma_n^2/M$，峰高近似源功率。
但峰附近二维偏差 $\boldsymbol\delta$ 下，若局部相关性非退化，写作

$$
\frac{|\mathbf a^H\mathbf a_0|^2}{M^2}\simeq1-\boldsymbol\delta^T\mathbf C\boldsymbol\delta,
\qquad P_C(\boldsymbol\delta)\simeq\frac{1}{p^{-1}+\boldsymbol\delta^T\mathbf B\boldsymbol\delta}
\quad(p\text{较大}).
$$

在固定局部二维区域积分，其主导质量随 $\log p$ 而非 $p$ 增长；一维非退化峰的对应局部质量则可随 $\sqrt p$ 增长。
这是有限孔径局部渐近分析，不是“所有二维算法必败”的证明。靠边界、角度退化、相关多源、具体加权时需另行分析。
光滑Jacobian可以改变局部权重，但不能普遍把这个峰面积自动校准为真实离散源功率。
还要注意，9阵元3×3半波长UPA的两轴孔径各为1个波长，而9阵元半波长ULA的线孔径为4个波长。
等阵元数并不意味着相同的方向分辨率；本轮未做等孔径控制，不能把现象纯粹归因于“二维”这个标签，更不能外推到所有大孔径UPA。

本轮固定第一邻区方向，只保留一个源、噪声方差1、无排除区、用总体协方差验证：

| 输入INR/dB | 真实功率p | 真方向Capon值 | dOmega谱的超白底积分质量 |
|---|---:|---:|---:|
| 0 | 1 | 1.1111 | 0.2169 |
| 10 | 10 | 10.1111 | 0.4970 |
| 20 | 100 | 100.1111 | 0.7931 |
| 30 | 1000 | 1000.1111 | 1.0838 |

超白底质量定义为 $\int(P_C-P_{\rm white})d\Omega=\operatorname{tr}(\mathbf R_{rec}-\mathbf R_{white})/M$，其中 $P_{\rm white}=1/M$。
不是把原始矩阵的整体尺度直接当成SINR误差：MVDR对整矩阵正比例缩放不敏感。
真正相关的是干扰质量相对白底的比例，以及连续展宽造成的协方差形状、特征子空间与零陷变化；单一标量缩放不能自动证明能够修正它们。

## 六、接收机白噪声不能直接当作各向同性远场散射

本模型的噪声是独立接收机噪声 $\sigma_n^2\mathbf I$。
若输入只有这个噪声，Capon谱是常数 $\sigma_n^2/M$，完整负半球立体角积分却得到

$$
[\mathbf R_{rec}]_{mn}=\frac{2\pi\sigma_n^2}{M}
\operatorname{sinc}\!\left(\frac{2\|\mathbf r_m-\mathbf r_n\|}{\lambda}\right).
$$

归一化sinc为 $\sin(\pi x)/(\pi x)$。
半波长UPA的对角邻居距离为 $\lambda/\sqrt2$，其归一化相关项约为 $\operatorname{sinc}(\sqrt2)=-0.217$，不是0。因此即便用整体尺度校准对角线，也不能变成白噪声单位阵。
扣掉球冠后还会进一步改变背景空间结构。

必须说明：半波长ULA在**各向同性三维散射的立体角模型**下有sinc整数零点；这不意味着原论文平坦 $d\theta$ 积分的白底也严格等于单位阵。两种测度不能偷换。

本轮用已知噪声白底矩阵作替换并按白输入迹校准，仅将SNR10 dB的立体角结果从9.4601提升到约9.7193 dB，仍远低于19.3608 dB。
这证明“把噪声底改对”不是已经验证的充分修复方案；连续谱展宽/功率质量映射仍需单独处理。
这里的白底替换使用已知噪声，是Oracle定位对照，不作为正式新算法。

## 七、权威资料支持什么，不支持什么

1. Björnson、Sanguinetti，IEEE Wireless Communications Letters，2021，*Rayleigh Fading Modeling and Channel Hardening for Reconfigurable Intelligent Surfaces*，DOI：10.1109/LWC.2020.3046107。
   [作者稿](https://arxiv.org/pdf/2009.04723)。正文第2页式(2)、(8)–(11)使用矩形平面阵列、含cos高度角的联合角密度和协方差二重积分，并得到sinc距离相关；第3页讨论矩形阵列对角方向的非零相关。
   它支持这里的方向测度和连续场相关模型，**不是原始Rec在当前场景一定高性能的证据**，也不是UAV INCM算法复现。

2. Nadeem等，*A Generalized Spatial Correlation Model for 3D MIMO Channels based on the Fourier Coefficients of Power Spectrums*，IEEE TSP作者稿，2015。
   [作者稿](https://arxiv.org/pdf/1504.05455)。第5页式(10)–(12)明确区分角度概率密度与含球面权重的功率谱定义，支持“密度必须和积分测度配套”。其阵列配置不是本工程的3×3 UPA，不能拿来冒充UPA Rec性能证明。

3. Gu、Goodman、Hong、Li，Signal Processing 96，375–381，2014，*Robust adaptive beamforming based on interference covariance matrix sparse reconstruction*，DOI：10.1016/j.sigpro.2013.10.009。
   [作者稿](https://arrc.ou.edu/~goodman/pubs/SP_robust_adaptive_beamforming.pdf)。第2–3页式(8)–(11)改用源支持选择与非负功率拟合，将接收机噪声单列为单位阵；第5页指出Capon重构的干扰功率低估问题。
   这是Gu后续工作的路线依据，不是当前CovFit与该论文逐步骤完全相同的证明，其ULA仿真损失数值也不能直接外推到本UPA。

本轮Tavily额度耗尽后改用普通检索读取作者公开正文，并目视核查平面阵论文第2页公式；不以摘要或标题推断积分定义。

## 八、复现入口与数据

```matlab
startup;
result=diagnose_upa_integral_mechanism("preliminary");
```

参数集中在 `config.m` 的 `cfg.experiments.upaIntegralMechanism`。
独立入口在 `experiments/diagnose_upa_integral_mechanism.m`，运行器为 `src/+uav/+experiment/runUpaIntegralMechanism.m`。
新增方位/高度角求积规则只扩展诊断核心 `src/+uav/+incm/integrateUpaBatch.m`，不切换正式方法。

结果前缀：`diagnosis_upa_integral_mechanism_20261108_205541`。
原始结果在 `results/raw/`，CSV在 `results/tables/`：`_summary`、`_power`、`_ablation`、`_contrast`、`_checks`。
最终检查与额外标度控制见下一节。数据带完整配置、代码版本、软件、时间与原始矩阵，时间明确为北京时间。

建议后续研究先定义期望恢复的功率测度及独立噪声模型，再研究密度校准/离散源重构与排除区设计；不能继续以“换一个Jacobian就一定正确”或“选一组参数让CovFit获胜”为目标。

## 九、最终验证与标度诊断

最终原始数据SHA256：`C90F3CFD755F64C59D761EE234404620AA16D2372EDF03FD008FF536CC5CC0FA`。
30项坐标、极细网格及5度排除区加密检查全部通过。
5度排除区3201→6401的五个SNR点中，立体角SINR最大变化0.0000735 dB；因此该对照的提升也不是偶然踩中某个网格节点。
单源质量3201→6401最大相对变化约3.75×10^-8，峰质量低估不能用漏采峰尖解释。

进一步只在SNR10 dB的总体控制组中测试

$$
\mathbf R_\beta=\beta(\mathbf R_{rec}-\mathbf R_{white})+\sigma_n^2\mathbf I,
\qquad \beta\in\{10^{-2},10^{-1},1,10,\ldots,10^6\}.
$$

在这些预设标度值中，du dv的最高SINR约12.0962 dB，立体角最高约10.4296 dB，均未达到19.3608 dB的真实上界。
所以当前问题不只是给干扰功率整体补一个常数：积分超白底矩阵的空间展宽/形状也会影响权值。
不能把这些离散标度扫描称为对所有实数标度的全局最优证明；它是对“仅补一个功率因子便可恢复”的反证线索。
这些标度只用于研究，不写入正式算法，不根据真值选择生产参数。

定向回归仅运行 `TestUpaIntegralBatch` 与 `TestUpaJacobianIntegral`，共5项测试全部通过。
MATLAB代码分析器对新运行器和修改后的积分核心无提示；`git diff --check`通过。
没有重跑实验6/7/8，没有修改期望导向估计或正式默认算法。
