function names = entryPoints()
%ENTRYPOINTS 返回五个公开实验入口的固定执行顺序。
%   run_all和测试共同调用此函数，避免各处复制不同顺序。

names = ["exp01_model_validation", "exp02_jitter_sweep", ...
    "exp03_snapshot_sweep", "exp04_dynamic_stability", "exp05_ablation"];
end
