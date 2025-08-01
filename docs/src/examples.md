
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
d = DistanceComoving(cosmology, 1000. * u"Mpc")
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

