# Experiment protocol

## Reproducibility and fairness

All runs start from `cfg=config(profile)`. Base seed `20260827` and a validated
substream scheme separate experiment, point, trial, scenario, period snapshot
and bootstrap streams. Each trial/period observation is generated once and
supplied unchanged to all methods or ablations. Method-specific seeds,
geometry, source powers or retuning are forbidden.

The fixed methods are Oracle, frozen, SCM, DL, static and temporal. Oracle alone
receives true INCM/ASV. Frozen retains first-period SCM-MVDR. SCM and DL use the
current SCM and nominal reference. Static and temporal use identical fixed
`gamma` and offline-calibrated `epsilonS`. Every same-period result checks the
Oracle SINR upper bound within numerical tolerance.

SINR is stored and aggregated in the linear domain, then converted to dB. The
95% interval is a deterministic paired cluster bootstrap: Monte Carlo trials
are resampled as the independent unit, all periods in a selected trial stay
together, and one trial-index draw is shared across methods.

## Profiles

- `smoke`: three trials, reduced sweep/ablation periods, full 20-period dynamic
  trace; pipeline validation only.
- `preliminary`: 100 paired trials and 100 periods/point.
- `paper`: 500 paired trials and 500 periods/point.

Main experiments use the configured 3x3 UPA; there is no array-size sweep.
Actual counts are stored in metadata.

## Experiments

1. Model validation exports attitude, every source direction cosine, normalized
   steering mismatch/array response, and a separate zero-jitter invariance test.
2. Jitter sweep reports median SINR with paired 95% intervals, 5th percentile,
   desired gain and interference leakage for all configured RMS points.
3. Snapshot sweep reports SINR, INCM NMSE, minimum eigenvalue, condition number
   and runtime for `L=[2,4,8,16,32,64]`.
4. Dynamic stability reports per-period SINR, empirical CDF, 5th percentile,
   outage probability, gain/leakage time series and mean/95th runtime.
5. Ablation compares full temporal reconstruction with one-at-a-time changes:
   `lambdaT=0`, identity alignment, uniform quadrature, no structure enhancement
   and pilot projection on. Full-method pilot projection remains off.

## Result schema and interpretation

Each experiment writes `<experiment>_<profile>.mat/.csv/.fig/.pdf/.png/.json`.
MAT contains calibrated cfg, code/environment metadata, substreams, raw
trial-period metrics, summaries, timestamp and runtime. PDF is vector; FIG is
editable; PNG is preview; CSV is compact summary; JSON is the manifest.
Overwrite occurs only when `cfg.output.overwrite=true`.

Do not use smoke plots to claim improvement. A defensible thesis statement
requires paired preliminary/paper intervals that support the difference across
the specified scans. Otherwise state only implementation completion and an
observed trend. Results do not generalize beyond this narrowband, fixed-source,
single-path cochannel model.
