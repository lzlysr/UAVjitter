# UAV attitude-jitter robust-INCM simulator

MATLAB R2024a simulation system for a UAV-mounted rectangular UPA receiving a
desired mmWave communication signal and fixed-world neighbouring-cell cochannel
interference. UAV attitude jitter changes the apparent array-frame directions;
output SINR is the primary metric.

The first version is narrowband, far-field and single-path. It deliberately
excludes 5G NR link waveforms, multipath, source motion, attitude-sensor input,
learning-based prediction and the final prediction/beamforming cascade.

## Requirements and run commands

- MATLAB R2024a or later;
- Optimization Toolbox for the robust spherical-mismatch solver;
- MATLAB Test for the automated suite.

From MATLAB:

```matlab
cd('D:\workspace\UAVjitter')
startup
r = runtests('tests', 'IncludeSubfolders', true);
assertSuccess(r)
run_all('smoke')
```

From a shell:

```powershell
matlab -batch "cd('D:\workspace\UAVjitter'); startup; r=runtests('tests','IncludeSubfolders',true); assertSuccess(r)"
matlab -batch "cd('D:\workspace\UAVjitter'); startup; run_all('smoke')"
```

`config(profile)` accepts `"smoke"`, `"preliminary"`, or `"paper"`.
`smoke` verifies the complete pipeline with reduced configured workloads;
`preliminary` uses 100 paired trials; `paper` uses 500. All physical,
algorithmic, workload, seed, plotting and export numerics are defined in
`config.m`.

## Implemented methods and experiments

The fixed comparison order is Oracle MVDR, frozen first-period MVDR, SCM-MVDR,
DL-MVDR, static Capon-INCM and temporal sparse structure-constrained INCM.
Static and temporal INCM use the same robust weight solver and one epsilon
calibrated offline from the configured design-jitter distribution. Non-Oracle
paths use only observations, configuration priors and previous estimated state.

The five entry points are `exp01_model_validation`, `exp02_jitter_sweep`,
`exp03_snapshot_sweep`, `exp04_dynamic_stability`, and `exp05_ablation`.
Each run writes a complete MAT result, CSV summary, editable FIG, vector PDF,
PNG preview and JSON manifest below `results/`. Generated artifacts and cache
files are ignored by Git.

## Evidence boundary

A passing `smoke` run establishes only that the model, algorithms, tests and
export pipeline execute consistently. It is not evidence that the proposed
method is superior. Statements such as “output SINR improved” require paired
`preliminary` or `paper` Monte Carlo results whose confidence intervals support
the difference. If they do not, report only implementation completion and an
observed trend. All conclusions remain conditional on this configured
narrowband cochannel-interference model and its parameter scans.

See [model_contract.md](docs/model_contract.md),
[variable_dictionary.md](docs/variable_dictionary.md), and
[experiment_protocol.md](docs/experiment_protocol.md).

## Layout

- `config.m`: sole production parameter entry;
- `src/+uav`: model, algorithms, metrics, experiment and export packages;
- `experiments`: five thin entry-point functions;
- `tests`: MATLAB Test contracts;
- `data/cache`: reproducible cache root;
- `results/raw|tables|figures|logs`: generated result bundle.
