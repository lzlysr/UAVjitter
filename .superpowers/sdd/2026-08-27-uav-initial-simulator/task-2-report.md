# Task 2 implementation report

## Scope

Implemented only Task 2 geometry, rectangular UPA, stationary attitude, and
snapshot-generation interfaces under `src/+uav`.  No beamforming, spectra, or
INCM functions were added.  The pre-existing tracked modification to
`task-1-report.md` was left untouched.

## RED evidence

Before production files existed, added `tests/TestTask2GeometryScenario.m` and
ran:

```powershell
matlab -batch "startup; r=runtests('tests/TestTask2GeometryScenario.m'); disp(r); assertSuccess(r)"
```

Result: `0 Passed, 8 Failed, 8 Incomplete` (exit status `0x00000001`).  Each
failure was the expected `MATLAB:undefinedVarOrClass` error for a missing Task
2 API: `rotationZYX`, `worldToArray`, `steeringUPA`, or `scenario.generate`.

## GREEN evidence

After minimal implementation, reran the same focused command.  Result:
`8 Passed, 0 Failed, 0 Incomplete` in 2.85 s.

Full current test suite command:

```powershell
matlab -batch "startup; r=runtests('tests','IncludeSubfolders',true); disp(r); assertSuccess(r)"
```

Result: `18 Passed, 0 Failed, 0 Incomplete` in 5.81 s.

## Delivered interfaces

- `uav.geometry.rotationZYX`, `worldToArray`, and `directionCoordinates`
- `uav.array.steeringUPA` with `vx kron vy` ordering for rectangular arrays
- `uav.scenario.generate` for stationary AR(1) attitude and steering truth
- `uav.scenario.generatePeriod` for seeded QPSK/noise snapshots, SCM, and
  algorithm-safe observations separated from `obs.truth`
- Focused coverage for rotation, zero attitude, rectangular steering, AR(1)
  covariance, determinism, `L < M`, `J = 0`, `J >= M`, finite dimensions, and
  large-snapshot SCM covariance convergence.

## Commit

The SHA of the commit that contains this report is recorded in the task
handoff/final delivery, because Git content addressing cannot embed that
commit's own final SHA in this file.

## Remaining concerns

The AR(1) covariance and SCM convergence checks are deliberately statistical
and use fixed streams with stated tolerances.  They verify the intended
distribution and covariance model without making later beamforming or INCM
design assumptions.

## Review-fix RED/GREEN evidence

The approved `obs = uav.scenario.generatePeriod(...)` interface remains a
single output: physical truth is intentionally nested under `obs.truth`.
This review fix does not replace it with a second output.

Added literal nonzero-ZYX, physical-input validation, QPSK independence, and
derived-noise-statistics checks, then ran:

```powershell
matlab -batch "startup; r=runtests('tests/TestTask2GeometryScenario.m'); disp(r); assertSuccess(r)"
```

RED result: `10 Passed, 1 Failed, 0 Incomplete` (exit status `0x00000001`).
The expected failing method was
`geometryAndSteeringRejectInvalidPhysicalInputs`: non-orthogonal rotations,
non-unit direction columns, and direction cosines outside the unit disk were
accepted without error.

After adding proper-rotation, unit-direction, and unit-disk validation, the
same command produced `11 Passed, 0 Failed, 0 Incomplete` in 3.74 s.  The new
signal test also verifies unit-modulus QPSK, low desired/interferer sample
correlation, and derived complex-noise power/circularity without exposing a
noise-truth field to algorithm inputs.
