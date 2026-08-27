# Task 1 configuration and execution contract

## Configuration

`cfg = config(profile)` accepts exactly `"smoke"`, `"preliminary"`, and
`"paper"`.  Every profile contains these sections:

`random`, `array`, `timeline`, `attitude`, `sources`, `signal`, `spectrum`,
`incm`, `beamforming`, `monteCarlo`, `experiments`, `plot`, and `output`.

`cfg.derived` is validation output, not an independent input source.  It holds
the element count `M`, interferer count `J`, signal/noise powers, world and
nominal array-frame unit vectors, nominal direction cosines, nominal steering
vectors, and the desired nominal steering reference.

The common physical defaults are:

| Item | Value |
| --- | --- |
| UPA | `Nx = Ny = 3`; spacing `[0.5, 0.5]` wavelengths |
| Desired direction | azimuth/elevation `[0, -35]` degrees |
| Interferers | `[40, -25]` and `[-55, -20]` degrees; INR `[20, 15]` dB |
| Desired-link SNR | `10` dB |
| Timeline | `dtSec = 1e-3`; `L = 8` snapshots per processing period |
| Attitude correlation | `[0.95, 0.95, 0.98]` for roll, pitch, yaw |
| Base seed | `20260827` |
| Jitter sweep | `[0, 0.5, 1, 2, 4, 8]` degrees |
| Snapshot sweep | `[2, 4, 8, 16, 32, 64]` |

Profiles differ only in workload/grid-resolution values: smoke uses 3 trials,
preliminary 100, and paper 500.  `K` and the lower-hemisphere spatial grid
increase with profile level.

Validation requires finite values, positive sizes/spacings/timestep/powers,
length-three attitude vectors, matching interferer direction/INR lengths, and
source elevations in the visible lower hemisphere `[-90, 0)`.  It supports no
interferers (`J = 0`) and arbitrary positive `J`.

## Reproducibility

`uav.util.makeStream(baseSeed, substream)` creates a local Threefry
`RandStream`; callers pass it explicitly to randomized production functions.
No production function may use or mutate global `rng` state.

## Startup and orchestration

`startup` adds only `src`, `experiments`, and `tests` under the repository root.
`run_all(profile)` executes the MATLAB Test suite first, then calls each
available `experiments/experiment_*.m` function in lexical order with the
validated configuration.  Its optional `RunTests` and `RunExperiments`
name-value switches exist only to make orchestration testable; both default to
true.

Generated cache and results are rooted at `cfg.output` and retained outside
version control using tracked `.gitkeep` placeholders.
