# UAV attitude-jitter INCM simulator

This repository contains a MATLAB R2024a simulator for a narrowband, far-field
UAV-mounted rectangular UPA receiving a desired mmWave communication signal and
fixed-world-frame neighbouring-cell cochannel interferers.  It evaluates output
SINR; it does not model 5G NR waveforms, source motion, multipath, learning, or
sensor input.

All numerical runtime parameters are defined in `config.m`.  Source directions
are supplied as azimuth/elevation degrees at that boundary; production geometry
uses radians.  Randomized production functions receive explicit `RandStream`
objects and do not alter MATLAB's global RNG state.

## Run

From MATLAB, open this repository and run:

```matlab
startup
r = runtests('tests', 'IncludeSubfolders', true);
assertSuccess(r)
run_all('smoke')
```

For a shell that starts MATLAB from another working directory, use the explicit
repository location:

```powershell
matlab -batch "cd('D:\workspace\UAVjitter'); startup; r=runtests('tests','IncludeSubfolders',true); assertSuccess(r)"
matlab -batch "cd('D:\workspace\UAVjitter'); startup; run_all('smoke')"
```

At the Task 1 stage, `run_all` runs the available test suite and discovers
`experiments/experiment_*.m` entry points in lexical order.  The five complete
experiment implementations and export products are intentionally deferred to
Task 5.

## Evidence boundary

The repository is a simulation framework, not experimental evidence for a
thesis claim.  Any output-SINR result is conditional on the configured
narrowband single-path cochannel-interference model, fixed world-frame sources,
and the selected profile.  Oracle quantities are reserved for future Oracle
methods and metrics, and must not enter non-Oracle algorithm inputs.

## Repository layout

- `config.m` — the sole source of numerical run and experiment parameters.
- `src/+uav` — package implementation, grouped by responsibility.
- `experiments` — future experiment entry points.
- `tests` — MATLAB Test contracts.
- `data/cache` and `results` — generated, git-ignored cache and output roots.

See [the Task 1 configuration contract](docs/contracts/task1-configuration.md)
for the current public fields and invariants.
