# Tavily 专项检索：IEEE 平面阵 INCM 重构的证据与复现边界

日期：2026-10-20。分支：`codex/gu-leshem2012-rec-study`。本报告补充而不覆盖此前本地 Zotero 核查记录。

## 1. 本轮实际做了什么

按照用户要求，实际启动本机配置的 `tavily-mcp@0.2.22`，通过标准 MCP stdio 客户端调用 `tavily_search` 与 `tavily_extract`，不是仅用其他搜索工具替代。最初免密连接发生限流；用户随后授权使用 API key，带 key 的连接已正常返回检索和抽取结果。

密钥通过隐藏输入进入临时子进程环境，不写入工程、配置或结果文件，不放在命令行参数中；结果保存前还做了凭据脱敏。检索完成后退出临时认证会话。无需为了本轮操作更改用户的 Codex/Tavily 配置。

本轮取得 **43 次搜索返回、2 次抽取返回**，共解析到 **438 条结果记录、316 个去重 URL**。这些数目包括误命中、同一论文的多个版本、非论文网页和不合格出处，**不是 316 篇论文，更不是 316 篇符合条件的论文**。另外保留了 13 次免密/接入排查阶段的限流记录，不计入上述返回数。

原始查询、结果与 URL 清单位于不提交 Git 的 `data/literature/tavily_upa_ieee_2026-10-20/`：

- `*_batch.json` / `*_search.json`：逐查询原始结果；
- `*_extract.json`：IEEE 官方页面和公开作者资料的抽取结果；
- `audit_summary.json`：查询状态与计数；
- `unique_urls.tsv`：网址去重清单，包含每条线索对应的查询和原始记录；
- `crossref_verified.json`：选定候选的官方题录交叉核查及错误记录。

Crossref 对 Lin 2021 返回 429，不能写成该条已获 Crossref 二次核验；其题录由本轮实际抽取到的 IEEE 官方页面确认。正文核查还使用了公开作者稿和作者单位页面；这部分访问与 Tavily 发现过程明确区分，不冒称全部由 Tavily 取得全文。

## 2. 严格筛选规则

优先关注 2019–2026 年，必要时追溯基础工作；搜索词覆盖 UPA、URA、uniform planar/rectangular、二维角域、Capon、INCM、covariance reconstruction/fitting、稀疏恢复、原子范数以及零陷展宽，并对关键题名和 IEEE 域名定向检索。

判断“符合 UPA 二维积分 INCM 证据”必须同时满足：

1. 正式 IEEE 期刊或领域会议，排除 IEEE Access；专业期刊、专门会议与旗舰期刊/会议分别注明，不一律当作同等级。
2. 正文确认物理阵列为均匀矩形平面阵，而非 ULA、圆阵、共形阵、稀疏平面阵或只有虚拟 UPA。
3. 目标是干扰加噪声协方差重构，不能只有 DOA 估计或普通接收协方差修复。
4. 取得重构公式，确认二维积分变量、积分区域和权重，并有对应数值验证。

仅标题有“planar”、在参考文献出现“URA”、二维 DOA 搜索、或只在引言提及 UPA，均不足以满足上述标准。Capon 谱估计中的二维扫描也不等于二维积分 INCM。

## 3. 检索结论先说明白

**本轮仍未确认满足全部四项要求的 IEEE-UPA 二维积分 INCM 论文全文。** 因此不能报告“已找到并复现某 IEEE-UPA 积分算法”，也不能反向声称“整个学界只做 ULA”或“二维积分在 UPA 下必然失效”。目前证据中，真实 URA 的非积分干扰重构、平面阵候选和稀疏极化阵功率拟合是三类不同材料。

与此前相比，本轮将最相关线索落到了可核对的官方题录、作者稿公式和实验设置，并明确了以下分级。

| 层级 | 论文 | 已核实事实 | 不能据此声称 |
|---|---|---|---|
| 真实 URA 非积分方法 | He 等，PAST 2022，*Two-Dimensional Adaptive Beamforming Based on Atomic-Norm Minimization*，DOI `10.1109/PAST49659.2022.9975052` | 作者稿为 6×6 URA；联合恢复干扰成分和目标方向，走原子范数优化，不做 Capon 角域积分 | 不是原始 Rec 的二维积分复现；PAST 是相控阵专门会议，不能写成 TSP/ICASSP 论文 |
| 最相关但未完成全文核查的候选 | Lin 等，IEEE Sensors Journal 21(4):5185–5194，2021，*Adaptive Beamforming Design of Planar Arrays Based on Bayesian Compressive Sensing*，DOI `10.1109/JSEN.2020.3030043` | 官方摘要确认平面阵、BCS 估计干扰 DOA、利用 DOA 和 Capon 谱重构，并加入旁瓣约束 | 物理阵是否严格均匀、重构是二维积分/局部积分/点源求和，均不能从摘要确定 |
| 功率拟合的直接 IEEE 方法论参照，不纳入严格 UPA 样本 | Yue 等，IEEE TVT 72(12):15648–15664，2023，*Adaptive Beamforming for Cascaded Sparse Diversely Polarized Planar Array*，DOI `10.1109/TVT.2023.3288611` | 作者稿式(46)–(50)给出协方差域最小二乘功率估计及干扰外积合成 | 物理阵为稀疏极化平面阵，含虚拟 CUPA 处理，不能冒充本工程普通 UPA；也不是当前 NNLS 的逐项复现 |

主要原始来源：[PAST 官方页面](https://ieeexplore.ieee.org/document/9975052)、[PAST 公开作者稿](https://www.researchgate.net/publication/362081194_Two-Dimensional_Adaptive_Beamforming_Based_on_Atomic-Norm_Minimization)、[作者单位发表记录](https://radio.seu.edu.cn/_s29/2018/1129/c19941a247869/page.psp)、[Lin 官方页面](https://ieeexplore.ieee.org/document/9220152)、[TVT 官方页面](https://ieeexplore.ieee.org/document/10159454)、[TVT 公开作者稿](https://www.researchgate.net/publication/371813879_Adaptive_Beamforming_for_Cascaded_Sparse_Diversely_Polarized_Planar_Array)。

## 4. PAST 2022：真实 URA，但不是积分路线

作者 Zeren He、Shengheng Liu、Xiaolong Miao、Yongming Huang；正式会议为 2022 IEEE International Symposium on Phased Array Systems & Technology，论文 1–5 页。作者稿第 II–III 节把二维阵列流形 Kronecker 向量化，随后从含期望信号的数据中分离目标项和稀疏干扰项，用原子范数及交替优化恢复干扰结构与目标方向。

第 IV 节使用 6×6 半波长 URA，目标二维角为 `(0°,0°)`，干扰为 `(0°,20°)`、`(0°,40°)`；波束图例使用 20 快拍、SNR 20 dB、每个 INR 10 dB；另做快拍/SNR 实验和 50 次蒙特卡洛。这与当前 3×3、8/64 快拍及三轴动态旋转不是同一场景，不能直接照搬性能结论。

对我们的启发是：二维阵列可以使用结构化稀疏恢复代替谱积分；但它不是“Capon 候选＋功率 NNLS”。若要做严格复现，还需核实二维 Toeplitz 变量的结构定义、方向更新的具体求解及干扰结构到可逆含噪协方差的处理。含未知方向的流形通常非线性，不能不作说明就把整个方向更新称为普通线性 SDP。

尤其是它会更新期望导向，而用户已经限定本工程先只研究 INCM。因此本轮不擅自把该联合估计器接入主实验；固定方向后仅借用其中一个环节，必须另标“组件改写”，不能继续称完整论文复现。文中的 `Traditional reconstruction` 曲线名称本身也不足以认定它就是本工程原始 Gu–Leshem Rec。

## 5. TVT 2023：与 CovFit 的联系落在功率拟合这一步

作者稿第 III-B 节式(46)的模型是：用估计的空间/极化导向组成矩阵 C，拟合

\[
\min_{\mathbf p}\left\|\widehat{\mathbf R}_x-\widehat\sigma_n^2\mathbf I
-\mathbf C\operatorname{diag}(\mathbf p)\mathbf C^H\right\|_F^2.
\]

向量化后，式(48)给出闭式最小二乘功率估计；式(50)重构时仅保留干扰对应项并加入噪声。也就是说，SCM 含期望信号并不禁止联合估计各分量功率，期望项在合成 INCM 时被排除。

当前 `upaCovfit.m` / `fitNonnegativePowers.m` 与此共享参数化协方差模型，但有三个实质差别：方向来自加载 Capon 候选而非其多参数求根；使用非负约束和可选 L2 正则而非闭式无约束最小二乘；物理阵、极化维度及可用先验不同。故可以引用它支持“协方差域估计功率后合成”的思路，不能写成当前代码完整复现该文。该文章也不能证明当前每一组 UPA 几何均有优势。

## 6. Lin 2021：摘要可读，不能据此补写未知公式

IEEE 官方页面明确显示 `You do not have access to this PDF`，正文在引言后要求登录。可确认 BCS 找干扰方向，再结合 Capon 谱重构；但未取得第 II–III 节和完整数值设置。

因此这里特别不做两个推断：既不把“使用 Capon”自动等同于全域二维积分，也不把“使用估计 DOA”自动等同于点源外积求和。这篇仍是最优先需要授权 PDF 的候选，当前没有合法全文，不能据摘要编写所谓原论文复现代码。

## 7. 其他重点线索的排除/降级记录

| 文献或线索 | 为什么不作为严格 UPA-INCM 积分证据 | 来源 |
|---|---|---|
| *Robust Covariance Matrix Estimation for Uniform Rectangular Array*，ICASSP 2026，DOI `10.1109/ICASSP55912.2026.11462357` | 确实是高质量会议及 URA；公开说明针对非高斯条件的块 Toeplitz 接收协方差估计，验证任务为 DOA，不能据此声称已经分离期望项重构 INCM | [会议官方程序](https://www.cmsworkshops.com/ICASSP2026/view_paper.php?PaperNum=14038&bare=1)；[作者仓储](https://centralesupelec.hal.science/hal-05663223v1/document)，本次直接访问被验证页拦截 |
| *2D Beamforming on Sparse Arrays with Sparse Bayesian Learning*，ICASSP 2019，DOI `10.1109/ICASSP.2019.8682747` | 官方摘要确认二维互质/嵌套稀疏阵和空间定位；不是物理均匀 UPA 上的 INCM-MVDR 重构验证 | [IEEE 官方页面](https://ieeexplore.ieee.org/document/8682747) |
| *Robust Wideband Beamforming of Planar Array for Rapidly Moving Interference*，ICSIDP 2019，DOI `10.1109/ICSIDP47821.2019.9173458` | 摘要确认平面阵及协方差锥削重构；这是宽带接收协方差与零陷展宽，不能等同于去期望分量的 INCM；会议层级也单独注明 | [IEEE 官方页面](https://ieeexplore.ieee.org/document/9173458) |
| *Linear Prediction-Based Covariance Matrix Reconstruction for Robust Adaptive Beamforming*，IEEE SPL 2021，DOI `10.1109/LSP.2021.3111582` | 官方摘要明确 ULA 双侧虚拟阵元扩展，虽为 IEEE 重构方法，仍不符合 UPA 要求 | [IEEE 官方页面](https://ieeexplore.ieee.org/document/9534656) |
| *Adaptive Interference Cancellation Using Atomic Norm Minimization and Denoising*，IEEE AWPL 2020 | Tavily 取得 NSF 公开稿；数值例明确为 4 阵元 ULA，不因“原子范数”便算作二维 UPA 证据 | [NSF 公开稿](https://par.nsf.gov/servlets/purl/10253411) |
| *Covariance Matrix Reconstruction-Based Multipath Suppression Beamforming for Uniform Rectangular Array*，ICSP 2026 | URA 相关，但属于多径、空间平滑/子空间路线；不是已确认的二维积分，不能以普通 IEEE 收录会议替代用户要求的核心会议层级 | [IEEE 页面](https://ieeexplore.ieee.org/document/11540684) |
| *Efficient covariance matrix reconstruction-based beamforming: a piecewise Clenshaw-Curtis quadrature approach* | 预印本存在二维 URA 求积描述，但没有本轮核实到的 IEEE 正式发表证据，不能纳入 IEEE 合格样本 | [作者预印本](https://www.researchsquare.com/article/rs-9730019/v1) |
| *Sparse Covariance Matrix Reconstruction-Based Nulling Broadening for UAV 2D Antenna Arrays* | 非 IEEE，且二维阵列不等于均匀矩形平面阵 | DOI `10.1155/2022/4987990` |
| *Two-dimensional adaptive beamforming for large planar array antennas based on weight matrix reconstruction* | IET RSN，而非 IEEE；题目中的 reconstruction 指权值矩阵，不自动属于 INCM 重构 | DOI `10.1049/rsn2.12151` |

另有大量 IEEE Access、普通阵列综合、医学成像、天线硬件设计及其他缩写误命中，全部不进入主证据列表。未取得正文的条目只按官方摘要作有限排除，不声称完整读过其算法。

## 8. 对当前 Rec、旧 static、CovFit 结论的修正

本轮没有新增能把三者差异归结为“UPA 二维积分天然错误”的证据。我们应继续区分：

- 原始 Rec：把 Capon 谱作为重构权重进行角域积分。
- 旧 static：在积分前后加入加载、谱底处理及结构修复，仍属积分路线。
- CovFit：选有限候选，在协方差域重新估计功率，再合成干扰项；不是只换积分网格。

下面是对现有代码模型的数学分析，而不是某篇论文的实验结论。令候选原子 \(\mathbf B_q=\mathbf a_q\mathbf a_q^H\)，其 Frobenius 内积为

\[
\langle\mathbf B_p,\mathbf B_q\rangle_F=|\mathbf a_p^H\mathbf a_q|^2.
\]

因此几何导致的流形相关性也会影响功率拟合的可辨识性和条件数。CovFit 并没有摆脱几何问题：候选漏检、相近候选、离网误差、噪声估计和少快拍的交叉项仍会影响结果。ULA 与 UPA 的阵列孔径、二维分辨能力及投影方向不同，也不能只根据一组曲线断言某种阵列天然更优。

最终可以写入研究记录的表述是：**已找到二维阵列非积分干扰恢复的 IEEE 先例，以及平面阵协方差域功率拟合的 IEEE 方法论依据；当前 Capon＋正则化 NNLS 仍是需独立验证的工程组合。** 这比“CovFit 无依据”和“CovFit 已被某篇普通 UPA 论文完整证明”都更准确。

## 9. 本轮是否实现新方法

没有。此前约定的条件是取得可核查的 IEEE-UPA 二维积分方法并据原文复现，目前尚未满足。不能把只读到摘要的 Lin 2021 改写成一个猜测实现；也不因找到 ANM 论文就越过“只研究 INCM，不加入期望导向估计”的范围。

本轮只增加检索工具和证据文档，没有修改 Rec、旧 static、CovFit、配置或开题报告，也没有启动新的 MATLAB 实验。新方法的源码和实验结果不存在，不能声称复现完成。

SI 状态：沿用用户“正文和 SI 都要”的选择，但本轮没有取得新的授权 PDF/SI 文件；作者稿网页可读与 PDF 已下载是两回事。无法访问的正文和补充材料仍标未取得，不写成“没有 SI”。

## 10. 复核与后续入口

本轮新增脚本均带中文说明：

- `scripts/tavily_literature_search.mjs`：实际 MCP 搜索/抽取与原始记录，遇服务端限流停止当前批次；
- `scripts/tavily_authenticated_session.py`：隐藏输入、临时环境认证，不持久保存 key；
- `scripts/summarize_tavily_literature.py`：汇总请求、网址去重和选定 DOI 的官方题录校验。

如继续推进，最优先的缺口仍是 Lin 2021 授权正文：核对物理阵均匀性、具体重构式、功率归一化及实验几何后，才能判断它是否真正触发原定复现流程。若只有抽象方法论可借鉴，则须明确标为本工程改写方案，使用配对几何/快拍扫描验证，不借用原论文整体算法名称或优越性结论。
