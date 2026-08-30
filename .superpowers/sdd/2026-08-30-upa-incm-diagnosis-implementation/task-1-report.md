# Task 1 实施报告：双快拍配置合同与三类UPA求积网格

## 完成范围

- 在 `cfg.literature.guLeshem2012.upaDiagnosis` 增加双快拍、四档网格、三类测度及诊断实验参数的唯一配置入口。
- 在 `validateGuLeshem2012` 中校验快拍、奇数网格、测度集合与选择项、加载、正则化、候选、抖动、试次、周期和实验编号。
- 新增 `uav.spectrum.upaIntegrationGrid`：方向余弦面积、方位角/z等立体角以及方位角/高度角立体角三类单元中心求积；流形只由配置与节点生成，不读取真值。
- 新增 Task 8 合同测试；按批准的最小范围更新既有 profile 等价性测试，清除新增的两个 profile 相关字段。

## RED

命令：

```text
matlab -batch "startup; r=runtests('tests/TestTask8UpaIncmDiagnosis.m'); assertSuccess(r)"
```

关键输出：两个测试均失败，错误标识为 `MATLAB:nonExistentField`，因为 `cfg.literature.guLeshem2012.upaDiagnosis` 尚不存在。失败位置分别为第 18 行和第 28 行；这证明新测试实际覆盖了缺失的配置合同。

## GREEN 与回归

命令：

```text
matlab -batch "startup; r=runtests('tests/TestTask8UpaIncmDiagnosis.m'); assertSuccess(r)"
matlab -batch "startup; r=runtests({'tests/TestTask8UpaIncmDiagnosis.m','tests/TestTask1Contracts.m'}); assertSuccess(r)"
```

关键输出：

- 聚焦测试：`2 Passed, 0 Failed, 0 Incomplete`。
- 最终聚焦加配置回归：`12 Passed, 0 Failed, 0 Incomplete`。
- `git diff --check` 无输出；三份新增或修改 MATLAB 文件经 `checkcode(..., '-id')` 无输出。

## 改动文件

- `config.m`
- `src/+uav/+config/validateGuLeshem2012.m`
- `src/+uav/+spectrum/upaIntegrationGrid.m`
- `tests/TestTask8UpaIncmDiagnosis.m`
- `tests/TestTask1Contracts.m`（经批准：归一化 `upaDiagnosis.monteCarloTrials` 与 `upaDiagnosis.periods`）

## 自审

- 未修改 `src/+uav/+incm/guLeshem2012Rec.m`。
- 新网格只读取 `cfg.array`、`cfg.spectrum` 和 `cfg.derived.nominalArrayUnitVectors`；源码检索未发现 `truth`。
- 三类网格的权重均为正且有限；`includedMask` 与 `desiredMask` 互斥，流形尺寸为 `M×N`。
- 配置中精确保留任务简报的指定数值，包括 `[8, 64]`、`[21, 51, 101, 201]`、三类测度和实验编号 `9`。

## 问题

无遗留问题。实施期间，既有 `TestTask1Contracts` 暴露新增的 `monteCarloTrials` 与 `periods` 都随 profile 改变；已按批准裁决将两者纳入其既有等价性归一化，未改变比较语义。

## 审查修复

审查指出方向余弦圆盘边界的 `du*dv/|uz|` 会令旧 `solid_angle` 测度振荡。已改为方位角/z参数化：z在[-1,0]取单元中心，`dOmega=daz*dz`，因此没有边界奇异，且总权重在浮点求和误差内等于 `2*pi`。`az_el_solid_angle` 保留为独立的方位角/高度角参数化。

同时将诊断合同收紧为：快拍数精确为`[8,64]`；加载为正有限标量；拟合正则化为非负有限标量；峰值阈值为`[0,1]`有限标量；NMS半径为非负整数标量；候选数、试次、周期、实验编号和选定网格均为正整数标量。选定测度继续要求为配置集合内的字符串标量。

修复RED命令：

```text
matlab -batch "startup; r=runtests('tests/TestTask8UpaIncmDiagnosis.m'); assertSuccess(r)"
```

修复前关键输出：旧 `solid_angle` 在21/51/101/201网格的总权重分别为 `6.08663658974816`、`7.65477121561286`、`6.29948604266582`、`6.66682815336498`，不等于 `2*pi`；零加载和多项标量字段未被拒绝，且零拟合正则化被错误拒绝。

修复GREEN与最终校验命令：

```text
matlab -batch "startup; r=runtests({'tests/TestTask8UpaIncmDiagnosis.m','tests/TestTask1Contracts.m'}); assertSuccess(r)"
matlab -batch "startup; checkcode('src/+uav/+spectrum/upaIntegrationGrid.m','-id'); checkcode('src/+uav/+config/validateGuLeshem2012.m','-id'); checkcode('tests/TestTask8UpaIncmDiagnosis.m','-id')"
git diff --check
```

最终关键输出：`14 Passed, 0 Failed, 0 Incomplete`；`checkcode`与`git diff --check`均无输出。21/51/101/201网格下，`az_el_solid_angle`总测度误差依次为`1.4650e-3`、`2.4836e-4`、`6.3324e-5`、`1.5989e-5`，单调下降；等立体角网格的排除后立体角误差分别为`4.7277e-2`、`3.1562e-3`、`4.6458e-4`、`1.5149e-3`，均受测试中的分辨率容差约束，并由单元中心对球冠边界的离散化误差解释。
