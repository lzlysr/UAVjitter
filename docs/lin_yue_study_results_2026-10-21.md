# Lin/Yue重构核心移植：快拍数实验与当前CovFit的性能差距

## 一、先给结论

本轮证据不支持“普通UPA下CovFit天然无效”，也不支持“正则5000就是唯一根因”。更准确的结论是：**当前实现把有限分辨率网格上的候选直接用于重构，候选的形状、位置精度与功率拟合相互耦合；增加快拍只能降低观测随机误差，不能自动消除这类字典误差。**

保持最终名义期望导向不变，使用完全相同的标准MVDR接口，仅以观测细化干扰候选后，非负拟合、Yue的无约束LS功率核心、Lin的Capon点功率核心在本场景都能恢复到接近19 dB的水平。它们在相同精细候选上十分接近，不能把这次恢复归功于某一个功率公式独有的优势。

这不是完整论文曲线复现：Lin的BCS及离网校正前端、旁瓣约束、期望导向更新未移植；Yue的不同极化稀疏平面阵和DSPR前端也未移植。此次是固定3×3普通UPA与固定最终期望导向条件下的**重构核心受控移植与根因诊断**。共同的亚网格Capon细化是明确新增的诊断前端，不冒称任何一篇论文的原算法。

## 二、结果来自哪里

- 分支：`codex/lin-yue-upa-snapshot-study`，基于main的`addef7b`；实现版本`7c68caf`。
- 工作树：`D:/workspace/UAVjitter-worktrees/lin-yue-upa-snapshot-study`，main和开题报告均未修改。
- 复放输入：`D:/workspace/UAVjitter/results/raw/exp07_classic_method_snapshot_sweep_preliminary.mat`，即用户指定FIG对应MAT，而非猜测当前config默认值。
- 输入实际参数：3×3半波长UPA，SNR=10 dB、两个INR=30/30 dB，期望方位/高度角[0,−35]°，干扰[40,−25]°、[−55,−20]°；三轴抖动RMS=5°，一阶相关系数[0.95,0.95,0.98]。
- 沿用归档20试次×20周期、51×51方向余弦网格；原扫描8:16:1024实际末点1016。本研究事前选取8个原扫描点：8、24、40、72、136、264、520、1016，未根据性能挑点。
- 另增加0°抖动配对诊断。相同原点序号、场景子流和快拍子流；不因方法改变随机实现。
- 原CovFit逐周期线性SINR与历史MAT的最大绝对误差为`2.44711e-11`，确认确实复现了用户所见问题。
- 所有INCM变体最终使用同一名义期望导向和`solve('mvdr',...)`；只有理想Oracle使用真实期望导向。DL/ESP/WCO仍保留各自经典基线定义，不把它们改造成标准MVDR。WCO仅引用相同5°观测的历史结果；0°诊断未计算WCO。

表中先汇总线性SINR的中位数，再转换dB。每个方法/快拍点最多400个周期样本，周期不是400次独立试验；区间按20个试次成簇bootstrap、每次500重采样。配对比较共用抽样，属于探索性逐对区间，未作多重比较校正。

## 三、5°目标场景：输出SINR随快拍数变化

单位：dB；完整8个快拍点见CSV。

|方法|L=8|L=72|L=1016|
|---|---:|---:|---:|
|理想Oracle|19.38|19.37|19.36|
|DL|5.51|13.41|17.33|
|ESP|5.51|14.56|18.71|
|WCO（同观测历史结果）|5.52|13.98|18.07|
|当前CovFit|13.04|15.50|15.44|
|同原候选，仅去L2正则|7.52|13.07|14.44|
|低加载局部峰＋零正则|10.79|11.19|10.93|
|亚网格干扰峰＋零正则NNLS|18.13|19.10|19.27|
|Yue-LS核心＋同一细化候选|18.13*|19.10|19.27|
|Lin原式核心＋同一细化候选|不适用|19.10|19.27|
|诊断：真实INCM＋固定名义导向|19.33|19.26|19.28|

`*`：Yue在L=8时3/400次重构INCM非PSD，此格只含397次成功样本；L=24时另有2/400次失败。不能把条件于成功的中位数当作无条件算法性能。NNLS细化版各点400/400均有有效结果并通过KKT审计。Lin在L=8<M=9时未加载SCM不可逆，原式标记不适用，没有暗中加载或用伪逆冒充原式。

相对当前CovFit，细化干扰峰＋NNLS的配对中位数增益：

- L=8：+5.09 dB，探索性95%区间[4.79,5.42] dB。
- L=1016：+3.83 dB，区间[3.45,4.10] dB。

这些结论仅适用于本轮几何、功率、网格和0°/5°抖动条件，不能外推为任意场景的优越性。

## 四、用受控实验逐项排查根因

### 4.1 不是这一场景下的最终期望导向失配主导

5°、L=1016时，真实INCM＋固定名义期望导向约19.28 dB，而原CovFit为15.44 dB；Oracle约19.36 dB。0°时原CovFit仍只有15.25 dB，细化版19.35 dB。因此本次约4 dB差距主要不是由最终期望参考偏差解释。**这不代表大抖动时不需要考虑期望失配；此前的大抖动性能边界仍然存在。**

### 4.2 最关键的对照：只细化干扰候选

低加载局部峰＋零正则与细化版共享加载0.001、期望字典、噪声估计、功率求解和最终MVDR，后者仅细化期望排除区外的干扰节点。5°、L=1016时由10.93提升到19.27 dB，配对增益8.33 dB，区间[7.87,8.78] dB。这是候选位置精度影响重构和零陷的直接证据，不涉及真实方向输入。

细化用固定邻域内的Capon分母最小化，不使用真实干扰数或真实方向。期望区候选不作为最终参考导向，也不执行期望导向校正。

设单位幅度阵列导向范数为sqrt(M)，候选与真实干扰的归一化相干度为

$$
\rho=\frac{|\widehat{\mathbf a}^{\mathrm H}\mathbf a|^2}{M^2}.
$$

将真实导向投影到候选的一维子空间，取最优复比例系数β，则

$$
\min_\beta\frac{\|\mathbf a-\beta\widehat{\mathbf a}\|_2^2}{M}=1-\rho.
$$

若权值仅在候选上形成零陷，真实干扰的剩余响应满足

$$
|\mathbf w^{\mathrm H}\mathbf a|^2
\leq M(1-\rho)\|\mathbf w\|_2^2.
$$

这只是误差上界，不是实际泄漏的等式。它说明“相干度接近1”不等于“足够支撑深零陷”。5°、L=1016时，每个样本对所有真实干扰的最小最佳相干度之中位数，原候选约0.997692，细化后约0.999999489。原候选在0.98覆盖阈值下已经覆盖两个干扰，**不能简单说它完全漏检了某个干扰**；问题主要是覆盖判据太宽松，位置误差对30 dB干扰仍有影响。

### 4.3 不能把问题简单归罪于正则5000

原候选＋原正则的干扰总拟合功率约为真实功率的5.3%；去正则后总功率约恢复到100%，但L=1016的SINR从15.44下降到14.44 dB，L=8下降更多。总功率正确不代表空间协方差形状正确。

诊断使用真实候选时，保留正则5000仍可达到约19.19 dB；去正则约19.28 dB。因此原正则会收缩功率，也会部分掩盖错误字典的影响，但不是本轮性能差距的唯一解释。这里的5000还依赖未归一化字典尺度；单个外积字典列的平方范数为M²=81，不能脱离流形归一化讨论其大小。

### 4.4 只降低加载、只做局部极大值筛选也不够

原加载0.12降至0.001、仍保留原正则时，L=8改善约1.44 dB，但L=1016反而下降约0.99 dB。只把相邻强网格单元合并成局部峰同样不能恢复性能。原因解释与实验相符：加载改变谱峰宽度和位置；固定网格局部峰依然存在偏差，相邻多个网格单元有时反而能近似一个离网方向。因此“峰更少、更像真实源个数”不能代替精度验证。

### 4.5 噪声估计不是本次主要瓶颈

5°、L=1016时，估计噪声功率中位数0.9757，真值1。仅替换成真实噪声，原CovFit约15.44 dB，没有恢复到19 dB。L=8<M时使用配置噪声功率回退，已明确记录，未冒称完全未知噪声估计。

## 五、论文公式和当前实现究竟哪里不同

### 5.1 Lin 2021：方向先定位，功率取Capon点值

正文5189页式（23）—（24）的核心为

$$
\widehat p_j=\frac{1}{\widehat{\mathbf a}_j^{\mathrm H}\widehat{\mathbf R}_x^{-1}\widehat{\mathbf a}_j},
\qquad
\widehat{\mathbf R}_{i+n}=\sum_{j\in\mathcal I}\widehat p_j\widehat{\mathbf a}_j\widehat{\mathbf a}_j^{\mathrm H}.
$$

这是有限候选求和，不是对二维Capon谱直接作全域积分。原式没有显式噪声单位阵项；本研究忠实保留这一点，最终MVDR只执行共同接口的最小数值加载。另设“加载Capon＋补噪声”的显式工程变体，但不称其为Lin原式。候选不准确时，直接在偏离真实峰的位置读取未加载Capon值可严重低估干扰功率；本轮相同原候选时L=1016只有约2.75 dB，换细化候选才恢复。

论文前端是BCS及离网校正，其20×20阵列具有与3×3完全不同的孔径和分辨率，还包含期望方向更新与旁瓣处理。不能把论文好结果仅归因于上面这一步功率求和。

### 5.2 Yue 2023：联合拟合，再剔除期望分量

正文15655页式（46）—（50）的普通阵列重构核心可写成

$$
\widehat{\mathbf p}=\arg\min_{\mathbf p\in\mathbb R^Q}
\left\|\widehat{\mathbf R}_x-\widehat\sigma_n^2\mathbf I
-\sum_{q=1}^{Q}p_q\widehat{\mathbf a}_q\widehat{\mathbf a}_q^{\mathrm H}\right\|_F^2,
\qquad
\widehat{\mathbf R}_{i+n}=\sum_{q\in\mathcal I}\widehat p_q\widehat{\mathbf a}_q\widehat{\mathbf a}_q^{\mathrm H}+\widehat\sigma_n^2\mathbf I.
$$

所有源（含期望）参与拟合，只有干扰项进入INCM。这里没有非负约束，也没有L2正则；代码用数值上更稳定的最小范数LS而不是显式求法方程逆。方向不是由LS自动分出“信号/干扰”，分类仍依赖固定期望先验区域和候选标签。

Yue使用不同极化稀疏平面阵、虚拟阵和DSPR前端。当前CovFit只继承了“从观测协方差联合估功率、剔除期望分量”的相似思想，并不等于Yue完整方法。用一堆相邻候选直接套无约束LS，会导致病态和负功率；本轮原候选LS产生非PSD矩阵4487/6400次，改用局部峰降为14次，细化候选降为6次。负功率但矩阵仍PSD的样本也没有暗中裁剪，其个数在原始结果中保留。

### 5.3 当前CovFit：在LS基础上增加非负和L2

$$
\widehat{\mathbf p}=\arg\min_{\mathbf p\geq0}
\left\|\widehat{\mathbf R}_x-\widehat\sigma_n^2\mathbf I
-\sum_q p_q\widehat{\mathbf a}_q\widehat{\mathbf a}_q^{\mathrm H}\right\|_F^2
+\lambda\|\mathbf p\|_2^2.
$$

这两项约束/正则是工程实现与Yue原式的差异；本研究没有为了使结果好看而隐去差异。同一精细候选上NNLS、Yue-LS和Lin点功率三者接近，说明当前首要研究对象应是**干扰候选误差与重构误差的耦合**，而不是仅继续盲调λ。

## 六、数值审计和不能隐藏的异常

对0°/5°、8个L、20×20全部样本，独立复放并核查非负二次规划KKT条件。令

$$
H_{qr}=|\widehat{\mathbf a}_q^{\mathrm H}\widehat{\mathbf a}_r|^2+\lambda\delta_{qr},\qquad
b_q=\Re\{\widehat{\mathbf a}_q^{\mathrm H}(\widehat{\mathbf R}_x-\widehat\sigma_n^2\mathbf I)\widehat{\mathbf a}_q\},\qquad
\mathbf g=\mathbf H\mathbf p-\mathbf b.
$$

检查p≥0、g≥0和p⊙g=0，使用归一化残差阈值1e−8。原CovFit与细化NNLS各6400个返回点全部通过，最大相对残差分别为2.30e−15、7.99e−16，排除这两条主要曲线由非负拟合未收敛导致的解释。

原候选零正则1/6400、低加载原候选零正则4/6400返回点不满足非负可行性，伴随`lsqnonneg`秩亏/迭代警告；这些返回点不能称为合格NNLS解。原始主结果未覆盖、未删样。单独报告剔除异常返回点的敏感性核查，受影响的4个组中位数最大改变0.0282 dB；这不是修复后的算法性能，也不足以解释本轮数dB差距。生产求解器的可行性/退出状态检查是后续应独立修复的问题。

无约束LS非PSD共4507次，均已保留异常信息；Lin在L=8不可逆、0°未计算WCO合计5600个“不适用/未计算”条目，不能混成求解失败率。

## 七、文件、参数与运行方式

全部新产物在此工作树`results/{raw,tables,figures,logs}/lin_yue_2026-10-21/`，不会覆盖main原图：

- 原始全方法MAT：`lin_yue_snapshot_study_preliminary.mat`，含实际cfg、元数据、子流、全部逐周期指标、失败状态和错误信息。
- 全方法汇总：`lin_yue_snapshot_study_preliminary.csv`。
- 配对效应与区间：`lin_yue_preliminary_paired_effects.csv`。
- 候选诊断：`lin_yue_preliminary_candidate_diagnostics.csv`。
- NNLS审计：`lin_yue_preliminary_nnls_audit_{detail,summary,sensitivity}.csv`及MAT。
- 主对比图：`lin_yue_preliminary_comparison.{fig,pdf,png}`。
- 0°/5°消融图：`lin_yue_preliminary_causal_ablation.{fig,pdf,png}`。
- 全方法四面板探索图：`lin_yue_snapshot_study_preliminary.{fig,pdf,png}`。主要结论以以上主对比和透明的失败率表为准。

本轮正式产物SHA256（重跑会覆盖研究目录同名文件，应先另行归档）：

```text
lin_yue_snapshot_study_preliminary.mat
2CC318CC18B0219F9FE76D4CA03A68B7C8F8F9AD76F275FDA08439D036D8B374
lin_yue_snapshot_study_preliminary.csv
95A58A5552659B23F4136DB0CE89D7A66B34BC2F335828B5DB0A532EAC10DA5C
```

在新工作树运行：

```matlab
startup;
% smoke：3试次×2周期；preliminary：归档20试次×20周期。
r = runLinYueSnapshotStudy("preliminary");
p = fullfile(pwd,"results","raw","lin_yue_2026-10-21", ...
    "lin_yue_snapshot_study_preliminary.mat");
analyzeLinYueSnapshotStudy(p);
auditLinYueNonnegativeFits(p);
analyzeLinYueAudit(p);
plotLinYueSnapshotStudy(p);
```

研究专用参数集中在`config.m`的`cfg.experiments.linYueStudy`。`pointIndices`是历史扫描点序号，不是快拍数本身；修改它时保持原序号才能准确复放。物理场景读取`sourceMat`归档cfg，不通过根config默认值暗改历史场景。

## 八、下一步边界

当前先不把细化版替换到main。建议独立验证更多几何、相近干扰、不同INR/快拍数和更大抖动，记录峰消失、合并及候选分类错误；再考虑BCS等论文前端的实质移植。本轮只验证一个目标几何与两个抖动条件，不能据此宣称完成两篇论文完整复现、性能普遍优于所有基线，或证明二维积分路线必然无效。

正文依据：Lin等，IEEE Sensors Journal 21(4),5185–5194，2021，DOI `10.1109/JSEN.2020.3030043`；Yue等，IEEE TVT 72(12),15648–15664，2023，DOI `10.1109/TVT.2023.3288611`。二者均已从本地Zotero正文核对；具体移植合同另见`lin_yue_study_protocol_2026-10-21.md`。
