# 本地 UAV jitter 文献库：UPA 与 INCM 重构路线核查

核查日期：2026-10-20。工作分支：`codex/gu-leshem2012-rec-study`。

## 1. 范围与证据边界

根据用户最新要求，本轮停止公网检索，只读取本机 Zotero 的 `UAV jitter` 分类（键 `TPNHRYML`）。本次本地 API 返回 54 个顶层条目，每个条目均有一个 PDF 附件，均取得全文索引，索引页数与总页数一致。这个计数不再沿用历史记录中的 22 篇。

完成的是全部条目的题录、摘要、全文关键词初筛，以及下述重点论文有关模型、重构与实验设置章节的核查；不代表逐字精读全部 54 篇。对公式缺失、双栏顺序错乱的全文索引，以本地 PDF 相关页核对，不能仅根据关键词或标题判定阵列类型。

本地清单与全文缓存：`data/literature/zotero_uav_jitter_2026-10-20/`。清单 `inventory.json` 保存题录、附件键、页数、全文长度和附件路径；缓存与受版权保护的 PDF 不提交 Git。本轮未修改 Zotero，不修改主算法，不启动 MATLAB 实验。

本分类返回的附件中未见独立 SI 附件；这仅表示当前分类附件列表没有，不表示出版商未提供 SI。本轮不联网补齐 SI。

## 2. 核心结论

1. 在本轮核查的本地库中，**尚未确认一篇同时满足“IEEE 期刊/会议、真实 UPA/URA、二维角域积分重构 INCM、给出对应数值验证”四项条件的论文**。这是本地样本的检索结果，不是对整个学界不存在该路线的证明，也不是二维积分必然不可行的证明。
2. 本库有明确的 UPA 姿态/信道建模论文，也有多篇明确的 ULA-INCM 论文；不能把二者拼起来说成已有 UPA 二维积分验证。
3. 方向/功率估计后重构协方差不是无来源的做法。本库存在 IEEE TAES 的参数化稀疏重构证据，也存在与当前 NNLS 功率拟合形式十分接近的非 IEEE 论文。**方法思想有依据，不等于当前 Capon 候选＋正则化 NNLS 实现已被某篇 IEEE-UPA 论文逐项验证。**
4. 不应把一篇联合导向估计、主瓣约束与 INCM 重构的整算法收益归因于重构本身。用户当前研究范围仍为 INCM 重构，不擅自引入期望导向校正。

## 3. 重点 INCM 论文分类

| 本地键 | 论文与出处 | 正文中的阵列/模型 | 重构路线与核查位置 | 对当前工程的意义 |
|---|---|---|---|---|
| `6XEHS7NA` | Gu–Leshem，2012，IEEE TSP；DOI `10.1109/TSP.2012.2194289` | 实验为 ULA | 排除期望扇区后积分 Capon 谱；另有期望导向估计 | 原始 Rec 基准，不是 UPA 二维验证 |
| `48PNC9FB` | Huang 等，2015，IEEE TSP；DOI `10.1109/TSP.2015.2396002` | 实验为 ULA | 对干扰导向矢量的环形不确定集合积分，离散求和逼近 | 积分域是导向不确定集合，不能误称 UPA 方位/俯仰二维积分 |
| `A7BTMI5R` | Zhang 等，2016，IEEE SPL；DOI `10.1109/LSP.2015.2504954` | 第 II 节 ULA | 第 III 节 SPSS：利用导向内积选择性、观测点采样与协方差锥削，避免密集 Capon 谱积分 | 可研究低成本采样，但依赖 ULA 结构和大阵列近似；不能不加检验迁到 3×3 |
| `FWYZTK2D` | Luo 等，2023，IEEE TAES；DOI `10.1109/TAES.2023.3263386` | 第 II 节 ULA；仿真 10 阵元，实测也是 ULA | 去除不需要的信号后，以 Gauss–Legendre 求积近似积分；另有导向估计 | 支撑非均匀求积的动机，不直接验证二维 UPA |
| `MGIK64S9` | Mohammadzadeh 等，2023，IEEE TAES；DOI `10.1109/TAES.2022.3232100` | 摘要及第 V 节明确 ULA，20 阵元 | PSEUR，第 III 节式(18)–(24)：两级分段常数谱＋干扰局部不确定扇区积分 | 并非“全域逐点硬积 Capon”；局部区域与功率模型同样重要 |
| `J6MR3LJB` | Mohammadzadeh 等，2023，IEEE TAES；DOI `10.1109/TAES.2023.3307507` | 单一空间频率 Vandermonde 流形；第 VIII 节 10 阵元 ULA | CMR-ISPS：最大熵谱、干扰邻域采样、加权外积重构、迭代权值计算 | 与局部谱采样有关，但不是现有 CovFit 的 NNLS，也不是 UPA 证据 |
| `9TJQUAKL` | Yue 等，2024，IEEE TAES；DOI `10.1109/TAES.2024.3397240` | 第 II 节 ULA，混合圆/非圆信号 | 第 III-C 节式(29)–(36)：增广协方差向量化、闭式最小二乘功率估计、只用干扰分量重构 | 明确支撑“方向与功率分开估计，再合成”的范式；不能照搬非圆信号增益到当前独立 QPSK |
| `2RI6E5WP` | Chen 等，2026，IEEE SPL；DOI `10.1109/LSP.2026.3692443` | 第 II 节 ULA，第 IV 节 16 阵元 | 式(7)、(10)–(13)：自定中心排除区＋加载 Capon＋减谱底＋一维积分；式(22)再加导数/WNG约束 | 与旧 static 最接近的新增参照，但仍不是 UPA 二维积分 |
| `JLQ5G8F6` | Chen 等，2026，Digital Signal Processing；DOI `10.1016/j.dsp.2026.106138` | 第 2.1 节及实验为 ULA | 式(25)–(28) SBL；式(33)–(43)门限、能量覆盖、方向/尺度解耦及外积合成 | 属于概率协方差模型与稀疏功率估计路线，不是 NNLS；非 IEEE，仅辅助 |
| `5PG9II8Y` | Chen 等，2026，Drones；DOI `10.3390/drones10030157` | 第 4.1 节，7 个 UAV 构成稀疏线阵，非 UPA | 第 3.3.2 节式(30)–(36)：联合期望/干扰 NNLS 功率拟合，重构时删除期望项 | 与当前功率拟合数学形式直接接近，但非 IEEE、依赖粗 DOA 与几何校准，不能替代高质量 UPA 实验证据 |

这些题录与正文来自用户本地附件，未通过公网重新核验元数据、出版状态或勘误。IEEE 的归类依据本地 PDF 页眉和 DOI，而不是联网验证结果。

## 4. Chen 2026 与原始 Rec、旧 static 的实际区别

原文 PDF 第 2–3 页（印刷页 2052–2053）给出：

\[
\widehat P_C(\theta)=\big[\mathbf a^H(\theta)(\widehat{\mathbf R}_x+c\mathbf I)^{-1}\mathbf a(\theta)\big]^{-1},
\]
\[
\widehat\theta\in\arg\max_{\theta\in\mathcal U(\theta_0)}\widehat P_C(\theta),\qquad
\theta_c=\begin{cases}\widehat\theta,&|\widehat\theta-\theta_0|>\tau,\\\theta_0,&\text{其他情况},\end{cases}
\]
\[
\widehat{\mathbf R}_{i+n}
=\int_{\Theta\setminus\mathcal U(\theta_c)}
\max(\widehat P_C(\theta)-\widehat\sigma_n^2,0)
\mathbf a(\theta)\mathbf a^H(\theta)\,d\theta
+\widehat\sigma_n^2\mathbf I.
\]

与原始 Rec 相比，新增加载、谱底扣除、自适应排除中心；与本工程旧 static 相比，同属“谱处理后积分”家族，但不能说实现相同：

- 旧 static 的期望排除区固定；原文排除中心可随局部 Capon 峰改变。
- 旧 static 在 UPA 方向余弦网格求积；原文用 ULA 单角度网格，步长 0.5°。
- 原文式(12)直接减去噪声功率；工程 `reconstruct.m` 使用 `(noiseEstimate + loading)/M` 作为白噪声谱底。单位幅度流形下，纯噪声的加载 Capon 谱为 `(σ²+c)/M`，所以两个式子不是同一个数值约定，不能不审查归一化就互换。
- 旧 static 最后做矩阵结构修复；原文的主要鲁棒性还包括导数约束和白噪声增益控制，不等价于普通 MVDR。
- 原文将期望角度建模为一维高斯偏差，干扰方向固定；不是本工程“三轴姿态旋转使所有世界固定来向共同改变”的动态模型。

因此，即便复现它在 ULA 的曲线，也不能据此解释本工程 UPA 的结果；整算法优于基线，也不能证明仅替换积分就足够。

## 5. 当前 CovFit 到底有哪些论文依据

### 5.1 IEEE TAES 2024：有明确的协方差域功率估计，但不是当前算法原样

Yue 等在第 III-C 节，将估计导向与非圆参数组成增广协方差字典，式(29)–(31)建立线性功率模型；式(32)采用

\[
\widehat{\mathbf p}=(\mathbf B^H\mathbf B)^{-1}\mathbf B^H\mathbf r.
\]

随后式(35)、(36)仅使用干扰对应的增广外积与噪声项合成 CAINC。这里的 \(\mathbf B\)、\(\mathbf r\) 位于共轭增广协方差域，不是当前工程的普通协方差字典。方向来自求根法，不是 Capon 峰；功率求解也没有当前 NNLS 的非负约束和 L2 正则项。

这是一条真实的 IEEE 方法论依据：先识别空间成分，再估功率，再只保留干扰成分。它不是“当前 Capon＋NNLS 在 UPA 下已有同一套实验证明”。

### 5.2 Drones 2026：与 NNLS 更直接接近，但严格降级为辅助材料

第 10 页式(30)–(35)明确写出

\[
\mathbf r=\operatorname{vec}(\widehat{\mathbf R}_x-\widehat\sigma_n^2\mathbf I),\quad
\mathbf G=[\operatorname{vec}(\widehat{\mathbf a}_s\widehat{\mathbf a}_s^H),
\operatorname{vec}(\widehat{\mathbf a}_{i,1}\widehat{\mathbf a}_{i,1}^H),\ldots],
\]
\[
\widehat{\mathbf p}=\arg\min_{\mathbf p\ge0}
\left\|\begin{bmatrix}\Re\mathbf G\\\Im\mathbf G\end{bmatrix}\mathbf p-
\begin{bmatrix}\Re\mathbf r\\\Im\mathbf r\end{bmatrix}\right\|_2^2,
\quad
\widehat{\mathbf R}_{i+n}=\widehat\sigma_n^2\mathbf I+
\sum_j\widehat p_{i,j}\widehat{\mathbf a}_{i,j}\widehat{\mathbf a}_{i,j}^H.
\]

它与本工程 `fitNonnegativePowers.m` 的基本目标一致；本工程改用保 Frobenius 范数的 Hermitian 实向量化，并可加 \(\lambda\|\mathbf p\|_2^2\)。但是其导向来自粗 DOA 先验与信标辅助几何校准，不是本工程二维 Capon 候选选择，且阵列不是 UPA。只能用作“NNLS 后剔除期望项”的直接形式参照，不能冒充满足用户筛选要求的 IEEE 论文。

### 5.3 为什么这仍不能消除 SCM 误差问题

拟合是在假设的物理模型集合上找近似协方差，不是在数学上自动消灭样本误差。若候选足够准确且数量受控，结构约束可能降低估计方差；若相邻候选高度相关、漏掉干扰、将期望泄漏放到排除区外，拟合也会分错功率。

期望和干扰的标签来自空间先验与候选区域，**不是 NNLS 自动认出哪一项属于用户信号**。联合拟合保留期望候选的作用，是给 SCM 中的期望能量一个解释位置；最后删除这些候选项才形成 INCM。若标签或候选有误，这种删除本身也会出错。残差小、KKT 条件满足或存在文献依据，都不保证 SINR 更优。

## 6. 本库明确涉及 UPA 的材料，应该用来支撑什么

| 本地键 | 论文 | 可以支撑 | 不能据此声称 |
|---|---|---|---|
| `ZJQ7EWNB` | Wang、Zhang，2022，IEEE TWC，`10.1109/TWC.2021.3118558` | UAV 姿态/位置到 UPA 二维角度及流形响应的几何映射 | UPA 下二维积分 INCM 已验证 |
| `J64533QB` | Zhao 等，2023，IEEE WCL，`10.1109/LWC.2023.3243590` | UAV 与地面 UPA、抖动进入阵列响应、ISAC 信道估计 | INCM 重构效果 |
| `5XXACXZ6` | Zhang 等，2019，IEEE ICC，`10.1109/ICC.2019.8761981` | UAV-to-UAV 的 UPA 坐标变换与姿态预测波束跟踪 | 纯阵列观测的无姿态输入 INCM 实验 |
| `K6V4SLI4` | Lee 等，2025，IEEE WCL，`10.1109/LWC.2024.3486785` | UPA 三维姿态与随机抖动下的波束优化 | 通过积分去除期望分量的重构方案 |
| `TR35PX5L` | 2025，IEEE CL，`10.1109/LCOMM.2025.3543284` | UPA 接收下角度不确定性、RSMA 鲁棒优化 | 由 SCM/Capon 重构 INCM |

其余阵列/算法相关条目的全文检索记录见本地清单。共形阵、多 UAV 分布式阵、稀疏线阵、只在引言提及 UPA 的论文均不自动纳入严格 UPA-INCM 证据。

## 7. 对后续实验的约束

当前未满足“找到 IEEE-UPA 二维积分方法并确认其验证”的条件，因此本轮不根据摘要或类比新增所谓论文复现方法，也不改变现有 Rec、旧 static、CovFit 定义。

后续仍可围绕已有代码研究：候选覆盖、离网误差、期望区分离、有限快拍功率估计误差与二维积分的模型/测度问题；但必须称为工程诊断或自主方案验证。不能预设“UPA 积分必然失败”或“CovFit 必然优于积分”。

如果只继续使用本地库，最有价值的下一步是分别建立两条证据链：UPA 姿态映射以 TWC 2022 为基础；INCM 估计以原始 Rec、URGLQ、PSEUR、参数化功率重构为参照。各自迁移到当前 3×3 UPA 时的新增假设和数学差异须单独记录，不将拼接后的方案命名为某篇原论文算法。

## 8. 复核方法

运行 `python scripts/audit_zotero_uav_jitter.py`，通过本地只读 API 重新生成分类清单与全文缓存。需要本地 Zotero 启动且允许本地 API；不需要 IEEE 网络权限。该脚本仅为资料清点，关键词命中不替代人工核查，不读取或写入 Zotero 数据库文件。
