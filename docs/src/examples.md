## Examples

Cosmos.jl ships with a handful of runnable scripts under the repository `examples/` directory. Each script demonstrates how one of the provided helpers behaves.

1. `examples/distance-conversions.jl` shows how to build `DistanceComoving`, convert it into `Redshift`/`ScaleFactor`, and move between distance measures without re-specifying the cosmology.
2. `examples/time-conversions.jl` walks through `TimeLookback`/`TimeConformal` and how to pass `Redshift`/`ScaleFactor` into the constructors.
3. `examples/custom-cosmology.jl` creates a non-Planck cosmology, highlights optional baryon handling, and prints the derived parameters.

To run an example:
```
julia examples/distance-conversions.jl
```

Each script is stand-alone and includes explanatory comments on the conversions being performed.
