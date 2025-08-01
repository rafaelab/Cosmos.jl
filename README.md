# Cosmos.jl

[![Build Status](https://github.com/rafaelab/Cosmos.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/rafaelab/Cosmos.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaelab.github.io/Cosmos.jl/index.html)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)


Useful tools for astrophysics and cosmology. 
It wraps and extends [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl/).

Things should work, but the code has not yet been thoroughly tested.


## Examples

### Planck 2018 cosmological parameters
```
using Cosmos

cosmology = CosmologyPlanck()
print(cosmology)
```

### Distance/Time Measures

It enables efficient conversion of distance and time measures in cosmology. 
```
d = DistanceComoving(cosmology, 600. * u"Mpc")
z1 = d |> Redshift 
z2 = convert(Redshift, d) # alternative
a = z1 |> ScaleFactor
D = convert(DistanceLuminosity, d)
```

It can also get a distance between two redshifts:
```
d = DistanceLightTravel(cosmology, Redshift(5.), Redshift(0.5))
```

Note that while these conversions might be convenient and intuitive, they are not (necessarily) efficient.


## Disclaimer

This program is provided 'as is', without warranties of any kind. 
Please use your discernment to interpret the results obtained with it.