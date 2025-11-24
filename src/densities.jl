# ----------------------------------------------------------------------------------------------- #
# 
@doc raw"""
	computeCriticalDensity(cosmology; z = 0.)

Computes the critical density of the universe at a given redshift:
```math
	\rho_\text{c} = 3 \dfrac{3 H^2(z)}{8 \pi G}
```

# Input 
- `cosmology::CosmologicalModel`: the cosmological model of interest 
- `z::Real`: the redshift at which to compute the density 
"""
@inline computeCriticalDensity(cosmology::CosmologicalModel; z::Real = 0.) = begin
	return upreferred(3. * hubbleParameter(cosmology, z) ^ 2 / 8π / G)
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc raw"""
	computeMatterDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
```math
	\rho_\text{m} = \Omega_\text{m} \rho_\text{c}
```

# Input 
- `cosmology::CosmologicalModel`: the cosmological model of interest 
- `z::Real`: the redshift at which to compute the density 
"""
function computeMatterDensity(cosmology::CosmologicalModel; z::Real = 0.)
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.Ωm * ρc)
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc raw"""
	computeRadiationDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
```math
	\rho_\text{r} = \Omega_\text{r} \rho_\text{c}
```


# Input 
- `cosmology::CosmologicalModel`: the cosmological model of interest
- `z::Real`: the redshift at which to compute the density 
"""
function computeRadiationDensity(cosmology::CosmologicalModel; z::Real = 0.)
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.Ωr * ρc)
end

# ----------------------------------------------------------------------------------------------- #
# 
@doc raw"""
	computeCurvatureDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
```math
	\rho_\text{k} = \Omega_\text{k} \rho_\text{c}
```


# Input 
- `cosmology::CosmologicalModel`: the cosmological model of interest
- `z::Real`: the redshift at which to compute the density
"""
function computeCurvatureDensity(cosmology::CosmologicalModel; z::Real = 0.)
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.Ωk * ρc)
end

# ----------------------------------------------------------------------------------------------- #
# 
@doc raw"""
	computeDarkEnergyDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
```math
	\rho_\text{\Lambda} = \Omega_\Lambda \rho_\text{c}
```


# Input 
. `cosmology::CosmologicalModel`: the cosmological model of interest
. `z::Real`: the redshift at which to compute the density
"""
function computeDarkEnergyDensity(cosmology::CosmologicalModel; z::Real = 0.)
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.ΩΛ * ρc)
end

# ----------------------------------------------------------------------------------------------- #
# 
@doc raw"""
	computeBaryonDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
```math
	\rho_\text{b} = \Omega_\text{b} \rho_\text{c}
```


# Input 
- `cosmology::CosmologicalModel`: the cosmological model of interest
- `z::Real`: the redshift at which to compute the density
"""
function computeBaryonDensity(cosmology::CosmologicalModel; z::Real = 0.)
	cosmology.Ωb < 0 || throw(ArgumentError("Cannot compute the baryon density because the baryon fraction was not provided to the `CosmologicalModel`."))
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.Ωb * ρc)
end


# ----------------------------------------------------------------------------------------------- #
# 