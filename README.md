# Cosmos.jl

[![Build Status](https://github.com/rafaelab/Cosmos.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/rafaelab/Cosmos.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaelab.github.io/Cosmos.jl/index.html)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Pragmatic helpers for cosmology that wrap [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl) and expose default conversions, rich distance/time value types, and typed helpers for common Planck-like models.

## Features

- `CosmologicalModel` keeps planck-compatible defaults, cached interpolation tables (distance/time ↔ redshift), and guards for optional parameters (e.g., baryons).
- Distance/time wrappers such as `DistanceComoving`, `DistanceLightTravel`, `TimeLookback`, and `TimeConformal` work with `Redshift`/`ScaleFactor` conversions out of the box.
- Typed conversion tables (`NamedTuple`s) keep lookups fast and avoid repeated `Dict` allocations.
- `CosmologyPlanck()` seeds the default cosmology, so you can build conversions without wiring up parameters manually.

## Getting started

```julia
using Cosmos

cosmo = CosmologyPlanck()
d = DistanceComoving(cosmo, 600u"Mpc")

redshift = d |> Redshift
scale = redshift |> ScaleFactor
luminosity = convert(DistanceLuminosity, d)
```

All distance/time conversions honour the cosmology stored in the value, so chained conversions stay consistent without re-specifying the model.

## Distance & time helpers

Along with the conversion wrappers shown above, you can create a distance from two redshifts or scale factors:

```julia
z1 = Redshift(2.0)
z0 = Redshift(0.1)
dt = DistanceLightTravel(cosmo, z1, z0)
```

This package pre-builds monotonic interpolations between redshift and any supported measure to keep repeated lookups fast, so repeated `DistanceXXX` constructors reuse cached tables rather than recomputing cosmology integrals.

## Tests & docs

- Run the bundled test suite with `julia --project -e 'using Pkg; Pkg.test()'`.
- Documentation is built with Documenter and lives at https://rafaelab.github.io/Cosmos.jl (see `docs/src` for the source).
- Example scripts live in `examples/`; each script is runnable via `julia examples/name.jl` to demonstrate common conversions.

## Contributing

Fixes and extensions are welcome. Please keep tests green, prefer descriptive doc strings, and run `Pkg.test()` before opening a PR.

## License

MIT © Rafael Abreu
