# 点Capon与局部高斯积分对照方法

## 定义及比较边界

2026-11-08（北京时间）在main新增两个独立可选方法，不替换现有方法。
当前UPA排除半宽5度是用户已有修改，保持不变；不为了新方法调整几何、功率或随机种子。
实验6、7、8默认组合追加两种方法，仍可自由选择任意子集；原主线实验1–5的方法合同不变。

两者的候选都来自当前观测Capon谱的局部极大值、相对门限、非极大值抑制及固定期望区排除。
默认候选参数与当前CovFit相同，但独立配置；峰数只受门限和最多M−1个候选限制，不读取真实J。
ULA共用 `cfg.spectrum.ulaGridSize`，UPA共用 `cfg.spectrum.gridSize`。
不作连续干扰方向细化，也不校正最终期望导向。

### incm_capon_point：离散源重构

$$
\widehat p_j=\frac{1}{\mathbf a^H(\widehat\Omega_j)\widehat{\mathbf R}_x^{-1}\mathbf a(\widehat\Omega_j)},
\qquad
\widehat{\mathbf R}_{i+n}=\sum_j\widehat p_j\mathbf a(\widehat\Omega_j)\mathbf a^H(\widehat\Omega_j)
+\widehat\sigma_n^2\mathbf I.
$$

不乘网格面积，不扣除峰点噪声底，不做NNLS。峰点功率仍可能受噪声、多源耦合、离网误差和有限快拍影响；不承诺接近上界。
噪声采用已有CovFit前端的小特征值中位数估计；L<M时按原合同回退到配置接收机噪声，不读取obs.truth。
前端为了接口复用会计算已有期望候选，但点重构只使用干扰部分，不拟合或调整期望方向。

### incm_urglq：局部角域积分与高斯求积

依据Luo等，IEEE TAES 59(5)，2023，DOI `10.1109/TAES.2023.3263386`，正文式(19)、(27)、(28)。
本工程只保留局部积分和求积环节：不构造期望抑制投影B，不采用投影后的准协方差，不做论文最后的期望导向校正。
直接用原SCM计算积分内的Capon值；因此必须称为“URGLQ重构环节对照”，不能称为完整论文性能复现。

ULA在候选等效入射角 `theta=asin(u_x)` 周围取正负8度的局部区间，裁剪到正负90度，积分测度为弧度dtheta。
默认三点Gauss–Legendre规则为

$$
\widehat{\mathbf R}_{i+n}\approx\frac{b-a}{2}
\left[\frac59\mathbf F(l_0)+\frac89\mathbf F(l_1)+\frac59\mathbf F(l_2)\right],
\quad
\mathbf F(\theta)=\frac{\mathbf a(\theta)\mathbf a^H(\theta)}
{\mathbf a^H(\theta)\widehat{\mathbf R}_x^{-1}\mathbf a(\theta)},
$$

其中 `l_n=(a+b)/2+(b-a)z_n/2`，标准节点为负sqrt(3/5)、0、正sqrt(3/5)。这些节点不是独立物理源。

UPA扩展已经用户确认：每个候选方向周围方位角、高度角各正负8度的矩形，每维3点，使用

$$
d\Omega=\cos\alpha\,d\phi\,d\alpha.
$$

下半球高度范围为负90至0度，方位角跨正负180度时周期拆分。
重叠ULA区间合并；重叠UPA矩形按边界分成不重叠单元，避免重复积分；因此总节点数不一定等于“峰数乘9”。
节点处继续施加固定名义期望排除区；当局部域碰到球冠边界时，低阶节点裁剪仅为近似，可提高求积阶数检查。
UPA扩展仍是同一Capon数值的立体角加权估计器，不代表已经证明它是原论文唯一正确推广。
高斯三点规则对低次多项式的代数精度，不等于对任意尖锐Capon峰都精确；需要参数扫描，不可凭名称宣称精度充分。

积分后默认不再显式加噪声I，以免悄悄改变积分定义；无候选或所有节点被期望排除时，采用估计噪声I回退，并记录fallback。
少量节点重构可能低秩，统一交给现有标准MVDR接口的最小数值加载处理，实际加载在beamforming诊断中记录。

## 配置与复现

唯一参数入口均在 `config.m`：

| 字段 | 默认值 | 作用 |
|---|---|---|
| `incm.caponCandidates.loadingRatio` | 0.001 | 只控制候选搜索谱加载 |
| `incm.caponCandidates.peakThresholdRatio` | 0.15 | 相对最强候选峰门限 |
| `incm.caponCandidates.peakNmsRadiusBins` | 0 | 局部峰间网格抑制半径 |
| `incm.caponCandidates.maximumCandidateCount` | 7 | 另受M−1限制的峰数上限 |
| `incm.caponPoint.powerLoadingRatio` | 0 | 点功率SCM加载比例 |
| `incm.urglq.powerLoadingRatio` | 0 | 积分内Capon的SCM加载比例 |
| `incm.urglq.quadratureOrder` | 3 | 每个角度维度的高斯节点数 |
| `incm.urglq.ulaHalfWidthDeg` | 8 | ULA局部角域半宽 |
| `incm.urglq.upaAzElHalfWidthDeg` | [8,8] | UPA局部方位/高度半宽 |
| `incm.urglq.upaMeasure` | solid_angle | 当前仅支持已确认的立体角扩展 |

两个功率加载比例为0时要求SCM可逆，L<M通常标记 `not-applicable`，MAT/CSV保留原因，曲线不伪造有效点。
候选谱的0.001加载不会隐式加到重构功率计算中；若确需运行秩亏场景，应明确将各自powerLoadingRatio设正数，并称为加载扩展。

```matlab
startup;
cfg=config("smoke");
names=["oracle","incm_integral_origin","incm_covfit", ...
    "incm_covfit_refined","incm_capon_point","incm_urglq"];
cfg.experiments.classicJitterMethodNames=names;
cfg.experiments.classicSnapshotMethodNames=names;
cfg.experiments.classicSnrMethodNames=names;
% cfg.array.Nx=9; cfg.array.Ny=1;  % 如需ULA，取消注释；否则默认3×3 UPA。
cfg=uav.config.validate(cfg);
r6=exp06_classic_method_jitter_sweep(cfg);
r7=exp07_classic_method_snapshot_sweep(cfg);
r8=exp08_classic_method_snr_sweep(cfg);
```

可按需只运行其中一行。默认结果仍写入 `results/raw`、`results/figures`、`results/tables`、`results/logs`；同profile/同阵列会沿用并覆盖相应实验文件名。
ULA输出附有 `_ULA_9x1` 等阵列标记。正式实验结果包含完整配置和原有版本/随机性元数据。

## 代码模块与验证

- `src/+uav/+incm/caponPoint.m`：离散源外积重构。
- `src/+uav/+incm/urglq.m`：局部积分重构，不进行额外结构修复。
- `src/+uav/+spectrum/interferenceCandidates.m`：共用候选与噪声前端。
- `src/+uav/+spectrum/caponAtDirections.m`：指定节点Capon计算与显式秩亏检查。
- `src/+uav/+spectrum/urglqQuadrature.m`：ULA/UPA局部域并集及高斯节点。
- `tests/TestCaponPointAndUrglq.m`：公式一致性、三点权重、球面面积、重叠去重、周期边界、真值隔离、秩亏行为、两个阵列及三个实验导出。

实验6/7/8始终使用相同名义期望导向和标准MVDR；不采用额外鲁棒接口，不进行方法专属真值调参。
本次仅进行定向测试和极小规模导出验证，不自动运行完整preliminary或替换既有正式结果。

本次新增7项测试均通过；`TestMatchedCovfit`与`TestAutomaticUla`回归通过，代码分析无提示。
另运行的 `TestIncmNamesAndSectors` 有3项旧测试未通过：它们硬编码UPA排除半宽25度，而用户当前config已改为5度。
未为通过这些断言回改用户参数，也未擅自改写旧测试；因此不能宣称整个仓库全部测试通过。
日志分别为 `results/logs/capon_point_urglq_tests_20261108.log`、`capon_point_urglq_regression_20261108.log` 与 `capon_point_urglq_codecheck_20261108.log`。
