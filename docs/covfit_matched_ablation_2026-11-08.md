# CovFit共用前端与纯细化消融

## 变更目的与基准

按用户要求，以现有`incm_covfit_refined`为唯一基准，不修改它的参数、期望字典、
干扰局部峰、连续细化目标及最终功率拟合定义。将`incm_covfit`替换为相同流程，
只关闭连续干扰方向优化。此变更同时适用于ULA与UPA，不修改积分型方法或temporal。

## 共用流程

相同观测与网格 → 相同低加载Capon谱 → 相同期望字典与干扰局部峰 →
可选连续干扰方向细化 → 相同非负功率拟合 → 相同KKT检查 → 剔除期望项 →
相同INCM合成 → 相同名义导向及标准MVDR。

“细化后处理”位于初始候选之后、最终功率拟合之前，并不是先拟合完功率再任意
移动方向。普通版完全跳过该连续优化；两种方法均不校正最终期望导向。
最终拟合功率与INCM允许不同，这是细化改变字典后的结果，不是更换了拟合方法或参数。

## 唯一参数入口

两种方法统一读取`cfg.incm.localCovfit`，保留原细化版数值：

| 字段 | 数值 | 使用范围 |
|---|---:|---|
| `loadingRatio` | 0.001 | 共用Capon前端及细化目标加载 |
| `fitRegularization` | 0 | 共用NNLS正则 |
| `peakThresholdRatio` | 0.15 | 共用候选门限 |
| `peakNmsRadiusBins` | 0 | 共用局部峰筛选 |
| `maximumCandidateCount` | 7 | 共用候选上限，另受M−1限制 |
| `kktTolerance` | 1e−8 | 共用功率解检查 |
| `radiusBins` | 2 | 仅连续细化分支 |
| `solver` | 沿用config.m | 仅连续细化分支 |

`cfg.incm.staticCovfit`保留为历史接口镜像；`validate`将五个共用字段从
`localCovfit`同步过去。请只通过`localCovfit`调整两方法的共同参数，不再用
`staticCovfit`为普通版单独设置加载或正则。

ULA四种静态重构仍共用`cfg.spectrum.ulaGridSize=2001`；UPA继续原profile网格。
所有网格、噪声及排除区规则在两种CovFit之间一致，不使用真实方向或真实干扰数。

## 文件与诊断

- `upaCovfit.m`：普通版公开入口，调用共用核心的关闭细化分支。
- `localCovfit.m`：共用完整实现，默认开启细化以保留原公开接口。
- `covfitFrontEnd.m`：共用期望候选与噪声估计，删除原来未参与结果的额外旧拟合。
- `refineInterfererCandidates.m`：只接收已选好的候选，不在细化分支重新选峰。
- `TestMatchedCovfit.m`：逐项核查两种方法的共享量及关闭细化分支的独立重构。

诊断保存`initialCandidateNodes`、`initialCandidateManifold`、`commonParameters`、
`noiseDiagnostics`和`refinementApplied`。同一观测下，前四项应逐元素一致。

## 历史结果与验证边界

旧普通版使用0.12加载、5000正则及最强网格单元候选；它与旧细化版的差异不只是
连续细化，因此旧图不能当作纯细化消融。当前普通版数值定义已改变，旧MAT/图
保留不改；若要复现旧普通版，应使用对应旧代码，而不是套用当前入口。

细化版的归档观测数值仍要求保持一致；归档测试只对细化版断言数值回归，普通版
仅记录相对旧基线的变化，不再错误要求它与旧值相等。此变更不预设谁必须获胜，
应在重新运行实验后同时查看中位数、低分位和失效比例。

本轮只做针对性测试及轻量导出，不运行正式preliminary，也不提交commit。

验证结果：42项针对性测试全部通过，包括两方法共享量逐元素相等、普通版独立NNLS
重构、UPA/ULA及8/32快拍、名称兼容、网格隔离及实验6/7/8轻量导出。原细化版
6组归档观测的线性SINR逐项一致；旧普通版差异属于本次主动更新，不再作为回归失败。
这不是全仓库测试通过声明，也不代表当前普通版或细化版在所有场景性能更优。
