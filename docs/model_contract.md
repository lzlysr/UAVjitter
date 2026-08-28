# Model contract

## Scenario boundary

The receiver is a UAV-mounted rectangular UPA in a 5G-cellular downlink
abstraction. One serving base station transmits the desired signal and `J`
normal neighbouring cells transmit cochannel interference. Source directions
are fixed in the world frame. Only UAV attitude changes their apparent
array-frame directions. This version has no source motion, barrage jammer,
multipath, fading, OFDM/DM-RS, coding, sensor measurement or learning model.

The equivalent complex-baseband model is narrowband, far-field and single-path.
Desired and interfering symbols are independent unit-power QPSK; noise is
circular spatially white complex Gaussian. Configured SNR/INR values set source
powers relative to `cfg.signal.noisePower`.

## Frames, angles and rotation

Configuration directions are `[azimuth,elevation]` in degrees. Elevation lies
in the visible lower hemisphere `[-90,0)` degrees. After validation, all angles
inside geometry and scenario functions are radians.

The attitude vector is `[roll,pitch,yaw]`. `rotationZYX` returns the active
array-to-world rotation `R = Rz(yaw) Ry(pitch) Rx(roll)`. A fixed world unit
direction is expressed in array coordinates as `u_array = R.' * u_world`. Its
direction-cosine coordinate is `xi=[u_array(1);u_array(2)]`. Proper-rotation,
unit-vector and visible-unit-disk invariants are validated.

## Attitude process and UPA

The independent attitude offsets follow
`DeltaEta_k=F_eta DeltaEta_(k-1)+nu_k`, with `F_eta=diag(rho)`,
`Q_eta=diag(stdRad.^2)` and `Q_nu=Q_eta-F_eta Q_eta F_eta.'`. The initial
offset is sampled from `N(0,Q_eta)`, so no burn-in is required. Every randomized
production function receives an explicit Threefry `RandStream`.

The `Nx` by `Ny` UPA lies in the array `x-y` plane. Steering entries have unit
magnitude and norm `sqrt(M)`, `M=Nx*Ny`. Sensor order is exactly `vx kron vy`:
the `y` index changes fastest within each `x` index. Rectangular `Nx~=Ny` is
supported.

## Observation and truth boundary

`generatePeriod` returns `obs.X [M x L]`, SCM `obs.Rx`, known desired pilot,
period index and nominal reference. Oracle/metric quantities are isolated under
`obs.truth`: true ASVs, true INCM, powers and directions. Non-Oracle
reconstruction is tested after removing or poisoning this field and may read
only observation fields, configuration priors and previous estimated state.

When `L<M`, physical noise is not identifiable from the SCM null spectrum in
the configured strong cochannel scene. INCM reconstruction therefore uses the
fixed configured receiver-noise prior and records
`configured-underdetermined`; for `L>=M` it uses a robust lower-eigenvalue
estimate. The Capon loading/noise pedestal is removed before source-spectrum
integration, then white noise is added once.

## Temporal INCM and robust weight

The temporal spectrum solves, without a hidden `1/2`,
`||p-c||^2+lambdaS||c||_1+lambdaT||c-GcPrev||_1`, `c>=0`. `G` is estimated by
bounded masked normalized 2-D correlation and becomes identity under weak or
ambiguous evidence. Peak count never uses true `J` and is capped at `M-1`.
Fine/coarse quadrature nodes remain in the unit disk and outside the nominal
desired exclusion.

The common robust weight minimizes `w^H RinHat w+gamma||w||^2` subject to
`real(w^H aRef)-epsilonS||w||>=1`, as a `2M`-real-variable Optimization Toolbox
problem. Epsilon is calibrated once from a fixed design-jitter distribution,
not from current trial truth. With `epsilonS=0`, it uses standard MVDR exactly,
independent of the configured positive `gamma` value.
