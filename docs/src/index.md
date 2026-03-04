# Cosmos.jl

Cosmos.jl wraps and extends [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl) with ergonomic value types, cached conversion helpers, and a Planck-motivated default cosmology.

Key additions:

- `CosmologicalModel` consolidates cosmology parameters, caching of conversion lookups, and optional baryon/neutrino metadata.
- Distance/time wrappers (`DistanceComoving`, `TimeLookback`, etc.) store their cosmology and expose friendly conversions toward `Redshift` and `ScaleFactor`.
- `CosmologyPlanck()` sets the default cosmology so you can skip parameter wiring for most use cases.
- Typed conversion tables (`NamedTuple`s) and `Unitful` helpers keep repeated lookups efficient and safe.

Refer to `docs/src/examples.md` or `examples/` to see runnable distance and time conversion snippets.
