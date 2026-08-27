# UAV Jitter Robust INCM Simulator Implementation Plan

> **For agentic workers:** Implement each task with strict red-green-refactor TDD. Do not dispatch subagents. Commit each completed task.

**Goal:** Build a reproducible MATLAB R2024a simulator for UAV attitude-jitter-induced array mismatch, cochannel interference, baseline beamforming, static Capon INCM reconstruction, and temporal sparse structure-constrained INCM reconstruction.

**Architecture:** `config.m` is the only numerical-parameter source. Pure package functions under `src/+uav` implement geometry, scenarios, algorithms, metrics, and plotting; experiment entry points only orchestrate configured sweeps. Truth is isolated from algorithm inputs and is available only to Oracle methods and metrics.

**Tech Stack:** MATLAB R2024a, MATLAB Test, Optimization Toolbox, optional Parallel Computing Toolbox. No CVX, no 5G waveform generation, and no learning subsystem.

**Spec:** User-approved plan in the conversation, grounded in `D:/workspace/hithesis-open/examples/hitart/reports/body/report_shenzhen_doctor_opening.tex` and `D:/workspace/hithesis-open/第三章研究内容及方案计划.md`.

## Global Constraints

- Work only in `D:/workspace/UAVjitter`; never modify the HiThesis report repository.
- Default main array is 3x3, while all production functions support rectangular `Nx ~= Ny` arrays.
- Interferer count `J` is derived from configured direction/power arrays; support `J=0` and arbitrary `J` without hard-coded loops or dimensions.
- The physical model is narrowband, far-field, single-path, fixed world-frame sources; only UAV attitude changes apparent array-frame directions.
- Current report scenario is normal neighboring-cell cochannel interference, not barrage jamming.
- All numerical runtime and experiment parameters live in `config.m`; tests may use explicit fixtures.
- Angles are degrees only at the configuration boundary and radians internally.
- Truth quantities never enter non-Oracle algorithm functions.
- MATLAB code uses explicit `RandStream` inputs for reproducibility and avoids global RNG state.
- No learning model, 5G NR waveform, multipath, source motion, sensor input, or report-source edits.

---

### Task 1: Repository Skeleton, Configuration, and Contracts

**Files:** Create `config.m`, `startup.m`, `run_all.m`, `README.md`, `src/+uav/+config/validate.m`, `src/+uav/+util/makeStream.m`, contract docs, output/cache directories, and focused tests.

**Requirements:**
- `config(profile)` accepts only `smoke`, `preliminary`, or `paper` and returns nested sections: random, array, timeline, attitude, sources, signal, spectrum, incm, beamforming, monteCarlo, experiments, plot, output.
- Defaults: Nx=Ny=3, half-wavelength spacings, desired az/el `[0,-35]` deg, interferers `[40,-25;-55,-20]` deg with INR `[20,15]` dB, SNR=10 dB, dt=1 ms, L=8, rho `[0.95,0.95,0.98]`, base seed 20260827, jitter sweep `[0,0.5,1,2,4,8]` deg, snapshot sweep `[2,4,8,16,32,64]`.
- Profiles vary only workloads/grid resolution: smoke 3 trials and short K, preliminary 100 trials, paper 500 trials.
- Validation derives M, J, powers, world unit vectors, nominal steering reference, and output paths; validates vector lengths, finite values, visible lower-hemisphere source directions, and positive sizes.
- `startup` adds only `src`, `experiments`, and `tests` paths.
- `run_all(profile)` runs tests then experiments in order.
- Tests must first fail because APIs are absent, then pass after implementation.

### Task 2: Geometry, Attitude, Rectangular UPA, and Snapshot Generation

**Files:** Create focused package functions under `+geometry`, `+array`, and `+scenario`, plus tests.

**Interfaces:**
- `R = uav.geometry.rotationZYX(etaRad)`.
- `uArray = uav.geometry.worldToArray(R,uWorld)` and `[xi,azElRad] = uav.geometry.directionCoordinates(uArray)`.
- `a = uav.array.steeringUPA(xi,Nx,Ny,dOverLambda)` with sensor ordering `vx kron vy`, unit-magnitude entries, and norm sqrt(M).
- `scenario = uav.scenario.generate(cfg,stream)` returns stationary AR(1) attitude offsets, rotations, directions, and steering truth without all snapshots.
- `obs = uav.scenario.generatePeriod(cfg,scenario,k,stream)` returns X, Rx, known QPSK pilot, algorithm-safe fields, and a separate truth structure containing true ASVs, true INCM, powers, and directions.

**Requirements:** stationary `Qnu=Qeta-F*Qeta*F'`; desired/interferer QPSK sequences independent; complex white noise; no fading; support rectangular arrays and J=0/J>=M. Test orthogonality, zero attitude, steering literals, stationary covariance, large-L covariance convergence, deterministic streams, L<M, and no NaN/Inf.

### Task 3: Baselines, Metrics, Capon Spectrum, and Static INCM

**Files:** Add focused functions under `+beamforming`, `+metrics`, `+spectrum`, and `+incm`, plus tests.

**Interfaces:**
- `[w,diag] = uav.beamforming.solve(method,R,aRef,cfg)` supports `mvdr` and `dlmvdr` here.
- `metrics = uav.metrics.evaluate(w,truth,cfg)` returns linear/dB output SINR, desired gain, interference leakage, and optional INCM NMSE.
- `[powerMap,diag] = uav.spectrum.loadedCapon(Rx,cfg)` scans the configured lower-hemisphere direction-cosine grid.
- `[Rhat,state,diag] = uav.incm.reconstruct(method,obs,priorState,cfg)` supports `static` here.

**Requirements:** implement Oracle MVDR, frozen initial MVDR orchestration support, SCM-MVDR, DL-MVDR, fixed nominal SOI exclusion, uniform quadrature, lower-eigenvalue noise estimate, Hermitian/PSD/shrinkage structure enforcement. Tests cover Capon peak localization, PSD/Hermitian output, Oracle SINR upper bound, J=0, and singular SCM loading.

### Task 4: Temporal Sparse INCM, Spectrum Alignment, Adaptive Quadrature, and Robust Weights

**Files:** Add temporal-spectrum, alignment, peak/refinement, structure-projection, and robust-weight functions with tests.

**Requirements:**
- Estimate a bounded integer 2-D spectrum shift using normalized correlation; use identity on weak evidence.
- Solve `||p-c||^2 + lambdaS||c||1 + lambdaT||c-GcPrev||1`, c>=0, by the exact separable piecewise proximal solution.
- Detect peaks without true J, cap at M-1, and construct coarse-outside/fine-inside non-overlapping quadrature nodes and area weights.
- Reconstruct, Hermitianize, eigenvalue-floor, and shrink the INCM.
- Optional pilot-orthogonal projection is an ablation switch and defaults off.
- Robust common INCM weight solves the real-form convex constraint `Re(w^H aRef)-epsilonS||w||>=1` with Optimization Toolbox; epsilonS=0 must return standard MVDR.
- Add `temporal` support to `uav.incm.reconstruct` and `robust` to beamforming solve.
- Tests cover known synthetic shifts, weak-evidence identity, proximal solution against hand-derived scalar cases and numeric search, non-overlapping quadrature, PSD output, robust feasibility, and epsilon-zero regression.

### Task 5: Experiments, Exports, Documentation, and End-to-End Verification

**Files:** Create five experiment functions, shared orchestration/aggregation/plot helpers, documentation, and integration tests.

**Requirements:**
- Exp 1: attitude/direction/steering-mismatch time series and zero-jitter validation.
- Exp 2: paired jitter sweep comparing Oracle, frozen, SCM, DL, static INCM, and temporal INCM; output median and 5th-percentile SINR, 95% confidence intervals, desired gain, and leakage.
- Exp 3: snapshot sweep with SINR, INCM NMSE, minimum eigenvalue, condition number, and runtime.
- Exp 4: dynamic SINR traces/CDF, low-percentile/outage, leakage, gain, and runtime.
- Exp 5: paired ablations for lambdaT=0, identity alignment, uniform quadrature, no structure enhancement, and optional pilot projection.
- Save complete cfg, MATLAB/toolbox metadata, method names, seeds/substreams, raw per-period/per-trial metrics, timestamps, and runtimes to MAT; export CSV summaries and FIG/vector-PDF/PNG figures.
- Smoke profile completes all five experiments and creates outputs. Generated outputs stay git-ignored.
- README documents the exact test and smoke commands and evidence boundaries for thesis claims.

**Final verification commands:**
`matlab -batch "startup; r=runtests('tests','IncludeSubfolders',true); assertSuccess(r)"`
`matlab -batch "startup; run_all('smoke')"`

