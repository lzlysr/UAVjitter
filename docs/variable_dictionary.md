# Variable dictionary

| Symbol / field | Size | Meaning | Access |
| --- | ---: | --- | --- |
| `Nx, Ny` | scalar | UPA element counts along array x/y | config |
| `M` | scalar | `Nx*Ny` sensors | derived config |
| `J` | scalar | interferer count inferred from direction/INR arrays | derived config |
| `K`, `L` | scalar | periods and snapshots per period | config / sweep |
| `eta`, `DeltaEta` | `3 x K` | `[roll;pitch;yaw]` attitude and offset, radians | scenario truth |
| `R_k` | `3 x 3` | array-to-world `Rz Ry Rx` rotation | scenario truth |
| `uWorld` | `3 x (J+1)` | fixed world-frame unit directions | config prior |
| `uArray_k` | `3 x (J+1)` | instantaneous array-frame directions | truth/diagnostics |
| `xi_k` | `2 x (J+1)` | array-frame x/y direction cosines | truth/diagnostics |
| `A_k` | `M x (J+1)` | true steering, desired first | Oracle/metrics |
| `aRef` | `M x 1` | nominal desired steering reference | algorithm prior |
| `X_k` | `M x L` | received complex snapshots | algorithm observation |
| `Rx_k` | `M x M` | sample covariance `X X^H/L` | algorithm observation |
| `Rin_k` | `M x M` | theoretical interference-plus-noise covariance | Oracle/metrics |
| `RinHat_k` | `M x M` | reconstructed INCM | algorithm estimate |
| `p_k`, `c_k` | grid matrix | residual Capon and enhanced temporal spectra | estimate/state |
| `G_k` | shift diagnostic | previous-to-current spectrum alignment | estimate/state |
| `w` | `M x 1` | complex beamforming weight | algorithm output |
| `sinrLinear` | scalar/array | desired output / interference-plus-noise | primary raw metric |
| `sinrDb` | scalar/array | `10*log10` after linear aggregation | reported metric |
| `desiredGain` | scalar/array | `|w^H a_s|^2` | metric |
| `interferenceLeakage` | scalar/array | total interferer output power, no receiver noise | metric |
| `incmNmse` | scalar/array | squared relative Frobenius error | metric |

Source columns are desired then configured interferers. Sensor order is
`vx kron vy`. Raw experiment arrays lead with `trial x period`, followed by
sweep point and method/variant where applicable. Snapshot signatures repeat
across methods/variants to prove paired observations. NaN is allowed only for
method-inapplicable INCM metrics; primary SINR/gain/leakage must be finite.
