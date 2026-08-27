# Task 1 report — Repository skeleton, configuration, and contracts

## Implementation summary

- Added the standalone MATLAB repository skeleton: `src/+uav`, `experiments`,
  `tests`, generated-cache/results roots, and tracked `.gitkeep` placeholders.
- Implemented `config(profile)` with the three approved profiles, immutable
  physical defaults, workload/grid-only profile changes, and a single numerical
  parameter source.
- Implemented `uav.config.validate(cfg)` to validate shape/finite/positive and
  lower-hemisphere constraints, derive array/source/power/steering quantities,
  and generate repository-relative output paths.  It allows zero or arbitrary
  numbers of interferers through derived `J`.
- Implemented explicit local Threefry stream creation in
  `uav.util.makeStream(baseSeed, substream)` without global RNG mutation.
- Added a narrow `startup` path initializer and `run_all` test-first,
  lexical-experiment orchestrator.  Per the Task 1 ruling, no experiment entry
  point exists yet; `run_all('smoke')` runs tests and discovers an empty set.
- Added README and a configuration/execution contract with the simulation and
  thesis-evidence boundaries.

## TDD record

### RED

Command run before the APIs existed:

```powershell
& 'D:\MATLAB\bin\matlab.exe' -batch "r=runtests('tests/TestTask1Contracts.m'); disp(r); assertSuccess(r)"
```

Expected result: exit code 1; `0 Passed, 7 Failed, 5 Incomplete`.
The failures were the expected absent contracts: `config` and `run_all` raised
`MATLAB:UndefinedFunction`, `uav.util.makeStream` raised
`MATLAB:undefinedVarOrClass`, and the absent `startup` path behavior did not
add the repository runtime paths.

### GREEN — focused contracts

```powershell
& 'D:\MATLAB\bin\matlab.exe' -batch "cd('D:\workspace\UAVjitter'); r=runtests('tests/TestTask1Contracts.m'); disp(r); assertSuccess(r)"
```

Result: exit code 0; `7 Passed, 0 Failed, 0 Incomplete`.

### GREEN — complete current suite

```powershell
& 'D:\MATLAB\bin\matlab.exe' -batch "cd('D:\workspace\UAVjitter'); startup; r=runtests('tests','IncludeSubfolders',true); assertSuccess(r)"
```

Result: exit code 0; `7 Passed, 0 Failed, 0 Incomplete`.

### Orchestrator smoke check

```powershell
& 'D:\MATLAB\bin\matlab.exe' -batch "cd('D:\workspace\UAVjitter'); startup; result=run_all('smoke'); assert(isempty(result.experiments))"
```

Result: exit code 0; the suite passed and the currently empty experiment set
was discovered without error.

## Files changed

- `config.m`
- `startup.m`
- `run_all.m`
- `src/+uav/+config/validate.m`
- `src/+uav/+util/makeStream.m`
- `tests/TestTask1Contracts.m`
- `README.md`
- `docs/contracts/task1-configuration.md`
- `data/cache/.gitkeep`
- `results/raw/.gitkeep`
- `results/figures/.gitkeep`
- `results/tables/.gitkeep`
- `results/logs/.gitkeep`

## Self-review

- Confirmed the approved physical defaults, profile trial counts, and profile
  invariance are exercised through behavior rather than source-text checks.
- Confirmed validation derives `M`, `J`, powers, world vectors, nominal steering
  reference, and output paths; it rejects invalid profile, source hemisphere,
  and direction/power length mismatch.
- Confirmed `makeStream` is deterministic for a seed/substream and leaves the
  global `rng` sequence unchanged.
- Confirmed `startup` has no recursive path addition and `run_all` executes
  tests before lexical experiment discovery.
- Ran `git diff --cached --check` before commit; it reported no whitespace
  errors.

## Concerns

- The local MATLAB installation runs a user-level startup that changes its
  initial directory.  README and verification commands therefore explicitly
  `cd('D:\workspace\UAVjitter')` before invoking repository scripts.
- The experiment directory is deliberately empty until Task 5.  This is the
  approved Task 1 boundary, not a failed smoke run.
