# Examples

Run any of the scripts via `julia examples/<name>.jl` to see the package helpers in action.

- `distance-conversions.jl`: shows typical distance conversions and how `Redshift`/`ScaleFactor` round-trips reuse a single cosmology.
- `time-conversions.jl`: demonstrates `TimeLookback`/`TimeConformal` conversions and unitful output.
- `custom-cosmology.jl`: builds a bespoke `CosmologicalModel`, prints it, and logs the associated matter and baryon densities.
