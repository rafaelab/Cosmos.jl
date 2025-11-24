# ----------------------------------------------------------------------------------------------- #
# 
@doc raw""" 
	scalingFunctionHubbleParameter(cosmology, z)

Compute the scaling function for the Hubble parameter (often dubbed `E`):  
```math
E(z) = \frac{H(z)}{H(0)} \sqrt{\Omega_r (1 + z)^4 + \Omega_m (1 + z)^3 + \Omega_k (1 + z)^2 + \Omega_\Lambda}
```
This follows the definition from Peebles 1993 (p. ~310-322), adopted by Hogg, arXiv:astro-ph/9905116.  


# Input
. `cosmo` [`CosmologicalModel`]: the cosmological model 
. `z` [`Real`]: the redshift 
"""
scalingFunctionHubbleParameter(cosmo::CosmologicalModel, z::Real) = Cosmology.E(cosmo.cosmology, z)
scalingFunctionHubbleParameter(cosmo::CosmologicalModel, z::Redshift) = scalingFunctionHubbleParameter(cosmo, z.value)
scalingFunctionHubbleParameter(cosmo::CosmologicalModel, a::ScaleFactor) = scalingFunctionHubbleParameter(cosmo, convert(Redshift, a))


# ----------------------------------------------------------------------------------------------- #
@doc """ 
	hubbleParameter(cosmology, z)

Compute the Hubble parameter H(z).  


# Input
. `cosmo` [`CosmologicalModel`]: the cosmological model
. `z` [`Real`]: the redshift
"""
hubbleParameter(cosmo::CosmologicalModel, z::Real) = Cosmology.H(cosmo.cosmology, z)
hubbleParameter(cosmo::CosmologicalModel, z::Redshift) = hubbleParameter(cosmo, z.value)
hubbleParameter(cosmo::CosmologicalModel, a::ScaleFactor) = hubbleParameter(cosmo, convert(Redshift, a))


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	hubbleConstant(cosmology, z)

Compute the present-day Hubble constant.


# Input
. `cosmo` [`CosmologicalModel`]: the cosmological model
"""
hubbleConstant(cosmo::CosmologicalModel) = hubbleParameter(cosmo, 0.)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	scaleFactor(z)

Compute the scale factor for a given cosmology. 
This should be the same for all cosmologies, by definition, but passing this argument fixes the correct type.


# Input
. `cosmo` [`CosmologicalModel`]: the cosmological model 
. `z` [`Real`]: the redshift 
"""
scaleFactor(cosmo::CosmologicalModel, z::Real) = eltype(cosmo)(1 / (1 + z))
scaleFactor(cosmo::CosmologicalModel, z::Redshift) = ScaleFactor(z)
scaleFactor(z::Real) = promote_type(Float64, typeof(z))(1. / (1. + z))
scaleFactor(z::Redshift) = ScaleFactor(z)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	hubbleDistance(cosmo)
	hubbleDistance(cosmo, z)

Computes the Hubble distance for a given cosmology and possibly at a given redshift.


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
"""
hubbleDistance(cosmo::CosmologicalModel) = Cosmology.hubble_dist0(cosmo.cosmology)
hubbleDistance(cosmo::CosmologicalModel, z::Real) = Cosmology.hubble_dist(cosmo.cosmology, z)
hubbleDistance(cosmo::CosmologicalModel, z::Redshift) = hubbleDistance(cosmo, z.value)
hubbleDistance(cosmo::CosmologicalModel, a::ScaleFactor) = hubbleDistance(cosmo, convert(Redshift, a))


# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	hubbleTime(cosmo)
	hubbleTime(cosmo, z)

Computes the Hubble distance for a given cosmology and possibly at a given redshift.  


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
"""
hubbleTime(cosmo::CosmologicalModel) = Cosmology.hubble_time0(cosmo.cosmology)
hubbleTime(cosmo::CosmologicalModel, z::Real) = Cosmology.hubble_time(cosmo.cosmology, z)
hubbleTime(cosmo::CosmologicalModel, z::Redshift) = hubbleTime(cosmo, z.value)
hubbleTime(cosmo::CosmologicalModel, a::ScaleFactor) = hubbleTime(cosmo, convert(Redshift, a))


# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	ageOfUniverse(cosmo)
	ageOfUniverse(cosmo, z)

Computes the Hubble distance for a given cosmology and possibly at a given redshift.  


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
"""
ageOfUniverse(cosmo::CosmologicalModel) = Cosmology.age(cosmo.cosmology, zero(eltype(cosmo)))
ageOfUniverse(cosmo::CosmologicalModel, z::Real) = Cosmology.age(cosmo.cosmology, z)
ageOfUniverse(cosmo::CosmologicalModel, z::Redshift) = ageOfUniverse(cosmo, z.value)
ageOfUniverse(cosmo::CosmologicalModel, a::ScaleFactor) = ageOfUniverse(cosmo, convert(Redshift, a))


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	comovingVolume(cosmology, z)

Calculates the comoving volume at a given redshift.  


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
"""
comovingVolume(cosmo::CosmologicalModel, z::Real) = Cosmology.comoving_volume(cosmo.cosmology, z)
comovingVolume(cosmo::CosmologicalModel, z::Redshift) = comovingVolume(cosmo, z.value)
comovingVolume(cosmo::CosmologicalModel, a::ScaleFactor) = comovingVolume(cosmo, convert(Redshift, a))


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	comovingVolumeElement(cosmology, z)

Calculates the comoving volume element at a given redshift.  


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
"""
comovingVolumeElement(cosmo::CosmologicalModel, z::Real) = Cosmology.comoving_volume_element(cosmo.cosmology, z)
comovingVolumeElement(cosmo::CosmologicalModel, z::Redshift) = comovingVolumeElement(cosmo, z.value)
comovingVolumeElement(cosmo::CosmologicalModel, a::ScaleFactor) = comovingVolumeElement(cosmo, convert(Redshift, a))


# ----------------------------------------------------------------------------------------------- #
# 
@doc raw"""
	comovingElement(cosmology, z)
	comovingElement(cosmology, a)

Calculates the comoving line element at a given redshift.
This is Hogg's eq. 28 adjusted.
```math
\frac{\mathrm{d} \ell}{\mathrm{d}z} = \frac{c}{H(z)} \frac{1}{E(z)}
```
NOTE: check nomenclature.  


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
- `z` [`Redshift`]: the redshift (of type `Redshift`) 
- `a` [`ScaleFactor`]: the scale factor type 
"""
comovingElement(cosmo::CosmologicalModel, z::Float64) = hubbleDistance(cosmo, z)  / Cosmology.E(cosmo.cosmology, z) / (1 + z)
comovingElement(cosmo::CosmologicalModel, z::Redshift) = comovingElement(cosmo, z.value)
comovingElement(cosmo::CosmologicalModel, a::ScaleFactor) = comovingElement(cosmo, convert(Redshift, a))


# # ----------------------------------------------------------------------------------------------- #
# # 
# @doc """
# 	calculateDensityNeutrinos(cosmology, Nν, Tcmb)

# Compute the density of relativistic components.

# # Input
# . `cosmo` [`CosmologicalModel`]: the cosmological model \\
# . `Nν`: effective number of neutrino species (3.04 according to Planck) \\
# . `Tcmb`: current CMB temperature \\
# """
# calculateDensityNeutrinos(cosmo::CosmologicalModel, Nν::Real, Tcmb::Real) = eltype(cosmo)((7. / 8.) * (4. / 11.) ^ (4 / 3) * Nν * calculateDensityPhotons(cosmo, Tcmb))
	

# # ----------------------------------------------------------------------------------------------- #
# # 
# @doc """
# 	calculateDensityPhotons(cosmology, Nν, Tcmb)

# Compute the density of photons.

# # Input
# . `cosmo` [`CosmologicalModel`]: the cosmological model \\
# . `Tcmb` [`Real` or `Temperature`]: current CMB temperature \\
# """
# calculateDensityPhotons(cosmo::CosmologicalModel, Tcmb::Real) = eltype(cosmo)(4.4813e-7 * (Tcmb ^ 2 / cosmo.h) ^ 2)


# ----------------------------------------------------------------------------------------------- #
# 