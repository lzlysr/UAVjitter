# 角域积分INCM重构是否只适用于ULA：检索与证据边界

检索日期：2026-10-20。本文件区分论文证据、代码实验和待验证推断，不以检索结果迎合预设结论。

## 1. 结论先行

“学界只在ULA上通过角域积分重构INCM，没有迁移到二维平面阵”不成立。已经核实二维平面稀疏子阵的Capon双重积分原式，也找到均匀矩形阵二维积分扩展的预印本。IEEE期刊中还存在面阵的稀疏参数估计、功率拟合和INCM重构路线，但不能把它们都叫作角域积分。

目前证据支持“二维重构的参数域、计算量和场景适配更复杂”，不支持“UPA性能必然不如ULA”。已有文献的成功场景和本项目的失败场景都不能直接证明一般性的阵型优劣。

## 2. 检索方式与限制

使用`nature-academic-search`的多源流程。当前会话未提供其专用学术MCP，使用技能的OpenAlex发现脚本、CrossRef官方DOI接口，以及网页检索、作者公开稿和出版社原文交叉核对。按DOI去重，优先核验原文的方法公式，而非只看标题中是否有“reconstruction”。未使用引文数或即时影响因子为论文排序，也不声称系统综述穷尽全部文献。

主要检索式：

- `"interference covariance matrix reconstruction" "planar array"`
- `"interference-plus-noise" "reconstruction" "uniform rectangular"`
- `"Capon" "covariance reconstruction" "two-dimensional"`
- `"planar array" "INCM" "integral"`
- `"covariance matrix reconstruction" "planar" "IEEE"`
- 对候选论文另作精确标题、DOI及作者稿检索。

部分IEEE页面不返回可读正文，MDPI/PMC出现429或验证页。相应条目明确标记摘要级或正文级，不把无法访问说成文献不存在。ResearchGate混有推荐文献，不能把页面下方其他论文的UPA描述归给当前论文。

## 3. 关键证据表

| 文献 | 出处与年份 | 已核验内容 | 对问题的作用 |
| --- | --- | --- | --- |
| Gu、Leshem：Robust Adaptive Beamforming Based on Interference Covariance Matrix Reconstruction and Steering Vector Estimation | IEEE TSP，60(7):3881–3885，2012 | 创始Capon角域积分路线；不是任意UPA性能定理 | 基础方法，不应由其ULA实验推导二维普遍成立 |
| Somasundaram、Jakobsson：Degradation of Covariance Reconstruction-Based Robust Adaptive Beamformers | SSPD，1–5，2014 | 作者原文：20阵元ULA中，干扰流形失配、离网格等也会使重构退化 | 说明失效机制并非UPA独有；其失配设置不同于本项目零失配几何实验 |
| Lin等：Adaptive Beamforming Design of Planar Arrays Based on Bayesian Compressive Sensing | IEEE Sensors Journal，21(4):5185–5194，2021；在线发表于2020 | DOI元数据及摘要级：面阵BCS估计干扰DOA，再结合Capon重构 | 面阵INCM确有研究，但不能据摘要将其归为二维连续积分 |
| Yang等：Robust Adaptive Beamforming Algorithm for Sparse Subarray Antenna Array Based on Hierarchical Weighting | Micromachines，13(6):859，2022 | 出版原文第7页式(28)：二维Capon双重积分；第8页式(33)还估计期望导向 | 直接否定“二维平面阵没有角域积分重构”；并非本项目3×3单层UPA的严格复现 |
| Yue等：Adaptive Beamforming for Cascaded Sparse Diversely Polarized Planar Array | IEEE TVT，72(12):15648–15664，2023 | 作者公开稿式(46)协方差Frobenius拟合，式(48)功率解，式(50)组成INCM | 高水平面阵重构证据；是稀疏点源参数路线，不是把Capon谱面积当功率 |
| Liu等：Efficient covariance matrix reconstruction-based beamforming: a piecewise Clenshaw-Curtis quadrature approach | Research Square，2026，v1预印本 | 第3.4节：12×12 URA，二维方向余弦、二维积分及张量积求积 | 明确的均匀矩形阵扩展线索，但未核实同行评审发表，不能包装为顶刊定论 |
| Ibrahim等：Robust Covariance Matrix Estimation for Uniform Rectangular Array | ICASSP，2026 | 官方摘要：非高斯环境下块Toeplitz协方差估计，用于DOA示例 | 高水平相关工作，但不是“排除期望区＋积分INCM”，故不计作该路线证据 |

### 3.1 可核查链接与DOI

1. Gu、Leshem，DOI：[10.1109/TSP.2012.2194289](https://doi.org/10.1109/TSP.2012.2194289)。[作者机构记录](https://cris.biu.ac.il/en/publications/robust-adaptive-beamforming-based-on-interference-covariance-matr-4/)。
2. Somasundaram、Jakobsson，DOI：[10.1109/SSPD.2014.6943304](https://doi.org/10.1109/SSPD.2014.6943304)。[作者公开全文](https://www.maths.lu.se/fileadmin/maths/personal_staff/Andreas_Jakobsson/SomasundaramJ14.pdf)。
3. Lin等，DOI：[10.1109/JSEN.2020.3030043](https://doi.org/10.1109/JSEN.2020.3030043)。[IEEE条目](https://ieeexplore.ieee.org/document/9220152/)。年份按正式卷期2021记录，不将DOI中的2020误当卷期年。
4. Yang等，DOI：[10.3390/mi13060859](https://doi.org/10.3390/mi13060859)。[出版社页面](https://www.mdpi.com/2072-666x/13/6/859)，[本次可读的出版PDF镜像](https://pdfs.semanticscholar.org/a610/f67fbb1330fdc7eb1ac609b75c0a34f07912.pdf)。
5. Yue等，DOI：[10.1109/TVT.2023.3288611](https://doi.org/10.1109/TVT.2023.3288611)。[IEEE条目](https://ieeexplore.ieee.org/document/10159454/)，[作者上传的公开稿](https://www.researchgate.net/publication/371813879_Adaptive_Beamforming_for_Cascaded_Sparse_Diversely_Polarized_Planar_Array)。
6. Liu等，DOI：[10.21203/rs.3.rs-9730019/v1](https://doi.org/10.21203/rs.3.rs-9730019/v1)。[预印本全文](https://assets-eu.researchsquare.com/files/rs-9730019/v1_covered_c079ab27-c25f-4f57-9183-792d7906cd1e.pdf)。
7. Ibrahim等，DOI：[10.1109/ICASSP55912.2026.11462357](https://doi.org/10.1109/ICASSP55912.2026.11462357)。[IEEE摘要](https://ieeexplore.ieee.org/abstract/document/11462357)，[会议官方记录](https://www.cmsworkshops.com/ICASSP2026/view_paper.php?PaperNum=14038&bare=1)。

## 4. 真正的二维积分证据与本项目差异

Yang等式(28)明确采用：

\[
\widehat R_{i+n}=\iint_{\overline\Theta}
\frac{a(\phi,\theta)a^H(\phi,\theta)}
{a^H(\phi,\theta)\widehat R^{-1}a(\phi,\theta)}\,d\phi\,d\theta .
\]

其引言认为当时一维应用多、二维较少，但论文自身已经做了二维扩展。“较少”不能改写为“没有”。其阵列为8个4×4平面子阵，先作子阵级固定波束形成，再在子阵间重构并估计期望导向；不是本项目的9阵元、固定名义期望参考。它采用平面角度测度，也不等同于本项目等立体角求积。存在这种论文既不能证明当前代码应自动成功，也不能被忽略后声称此路线从未迁移。

Liu等预印本第3.4节进一步明确写出URA的`a_x(u)⊗a_y(v)`，并讨论二维积分、块Toeplitz结构和张量积求积。它主要提供大阵列方向图扩展，不能当成当前小UPA几何扫描的严格对照，更不能因网页称其为“Research Article”就断言已同行评审。

## 5. CovFit的相关论文依据与不可混同之处

Yue等式(46)在估计各源参数后拟合：

\[
\min_p\|\widehat R_{xx}-\widehat\sigma_n^2I-C\operatorname{diag}(p)C^H\|_F^2,
\]

并用式(50)只保留干扰项组成INCM。这与本项目“点源字典＋协方差拟合＋去掉期望项”的建模思路相关，因此不能说协方差拟合毫无文献基础。但其方向/极化估计、阵列结构及无约束闭式功率解均与我们的Capon候选、非负最小二乘、可选正则和局部细化不同。**本次找到的是相关依据，不是把当前代码追认成该论文复现。**

这篇TVT还比较多维积分方法并指出其计算负担；计算效率问题与“UPA下SINR必然更低”是两件事。本文未引入它的期望导向估计，仍遵守用户只研究INCM的边界。

## 6. 新发现的高度相关论文：需另行精读，不冒用为UPA证据

Chen等，Robust UAV Beamforming Under Pointing Jitter Using Self-Centered Covariance Reconstruction，IEEE Signal Processing Letters，33:2051–2055，2026，DOI：[10.1109/LSP.2026.3692443](https://doi.org/10.1109/LSP.2026.3692443)，[IEEE条目](https://ieeexplore.ieee.org/document/11516153/)。

本次已通过CrossRef核实题名、卷页、作者和DOI。其主题与本项目高度接近，值得优先获取并精读；但未取得可核验的全文阵型和公式，不能凭检索页面混入的其他论文摘要，断言它采用UPA或将其算作二维积分的直接证据。未联系作者、未修改Zotero或开题报告。

## 7. 当前适合写入研究记录的表述

> 既有INCM角域积分方法大量围绕一维阵列展开，已有研究将重构扩展到二维平面阵或采用面阵稀疏参数重构。对于本项目的小规模UPA、有限快拍和姿态抖动场景，直接二维化的重构性能表现出明显的几何、测度及角域分区敏感性，因此有必要开展针对性诊断和改进。

不宜写：“二维UPA没人做”“原论文只对ULA有效”“角域积分在UPA下必然失败”或“文献数量能证明UPA性能低于ULA”。独立ULA实验见[ULA几何诊断报告](rec_ula_geometry_study_2026-10-20.md)。
