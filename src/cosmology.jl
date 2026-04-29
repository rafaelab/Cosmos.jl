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
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
"""
function scalingFunctionHubbleParameter(cosmo::CosmologicalModel, z::Real)
	return Cosmology.E(cosmo.cosmology, z)
end

function scalingFunctionHubbleParameter(cosmo::CosmologicalModel, z::Redshift)
	return scalingFunctionHubbleParameter(cosmo, z.value)
end

function scalingFunctionHubbleParameter(cosmo::CosmologicalModel, a::ScaleFactor)
	return scalingFunctionHubbleParameter(cosmo, convert(Redshift, a))
end


# ----------------------------------------------------------------------------------------------- #
#
@doc """ 
	hubbleParameter(cosmology, z)

Compute the Hubble parameter H(z).  


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model
- `z` [`Real`]: the redshift
"""
function hubbleParameter(cosmo::CosmologicalModel, z::Real)
	return Cosmology.H(cosmo.cosmology, z)
end

function hubbleParameter(cosmo::CosmologicalModel, z::Redshift)
	return hubbleParameter(cosmo, z.value)
end

function hubbleParameter(cosmo::CosmologicalModel, a::ScaleFactor)
	return hubbleParameter(cosmo, convert(Redshift, a))
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	hubbleConstant(cosmology, z)

Compute the present-day Hubble constant.


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model
"""
function hubbleConstant(cosmo::CosmologicalModel)
	return hubbleParameter(cosmo, 0.)
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	scaleFactor(z)

Compute the scale factor for a given cosmology. 
This should be the same for all cosmologies, by definition, but passing this argument fixes the correct type.


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
"""
function scaleFactor(cosmo::CosmologicalModel, z::Real)
	return eltype(cosmo)(1 / (1 + z))
end

function scaleFactor(cosmo::CosmologicalModel, z::Redshift)
	return ScaleFactor(z)
end

function scaleFactor(z::Real)
	return promote_type(Float64, typeof(z))(1. / (1. + z))
end

function scaleFactor(z::Redshift)
	return ScaleFactor(z)
end


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
function hubbleDistance(cosmo::CosmologicalModel)
	return Cosmology.hubble_dist0(cosmo.cosmology)
end

function hubbleDistance(cosmo::CosmologicalModel, z::Real)
	return Cosmology.hubble_dist(cosmo.cosmology, z)
end

function hubbleDistance(cosmo::CosmologicalModel, z::Redshift)
	return hubbleDistance(cosmo, z.value)
end

function hubbleDistance(cosmo::CosmologicalModel, a::ScaleFactor)
	return hubbleDistance(cosmo, convert(Redshift, a))
end


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
function hubbleTime(cosmo::CosmologicalModel)
	return Cosmology.hubble_time0(cosmo.cosmology)
end
function hubbleTime(cosmo::CosmologicalModel, z::Real)
	return Cosmology.hubble_time(cosmo.cosmology, z)
end
function hubbleTime(cosmo::CosmologicalModel, z::Redshift)
	return hubbleTime(cosmo, z.value)
end
function hubbleTime(cosmo::CosmologicalModel, a::ScaleFactor)
	return hubbleTime(cosmo, convert(Redshift, a))
end


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
function ageOfUniverse(cosmo::CosmologicalModel)
	return Cosmology.age(cosmo.cosmology, zero(eltype(cosmo)))
end
function ageOfUniverse(cosmo::CosmologicalModel, z::Real)
	return Cosmology.age(cosmo.cosmology, z)
end
function ageOfUniverse(cosmo::CosmologicalModel, z::Redshift)
	return ageOfUniverse(cosmo, z.value)
end
function ageOfUniverse(cosmo::CosmologicalModel, a::ScaleFactor)
	return ageOfUniverse(cosmo, convert(Redshift, a))
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	comovingVolume(cosmology, z)

Calculates the comoving volume at a given redshift.  


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
"""
function comovingVolume(cosmo::CosmologicalModel, z::Real)
	return Cosmology.comoving_volume(cosmo.cosmology, z)
end
function comovingVolume(cosmo::CosmologicalModel, z::Redshift)
	return comovingVolume(cosmo, z.value)
end
function comovingVolume(cosmo::CosmologicalModel, a::ScaleFactor)
	return comovingVolume(cosmo, convert(Redshift, a))
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	comovingVolumeElement(cosmology, z)

Calculates the comoving volume element at a given redshift.  


# Input
- `cosmo` [`CosmologicalModel`]: the cosmological model 
- `z` [`Real`]: the redshift 
"""
function comovingVolumeElement(cosmo::CosmologicalModel, z::Real)
	return Cosmology.comoving_volume_element(cosmo.cosmology, z)
end
function comovingVolumeElement(cosmo::CosmologicalModel, z::Redshift)
	return comovingVolumeElement(cosmo, z.value)
end
function comovingVolumeElement(cosmo::CosmologicalModel, a::ScaleFactor)
	return comovingVolumeElement(cosmo, convert(Redshift, a))
end


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
function comovingElement(cosmo::CosmologicalModel, z::Real)
	return hubbleDistance(cosmo, z) / Cosmology.E(cosmo.cosmology, z) / (1 + z)
end
function comovingElement(cosmo::CosmologicalModel, z::Redshift)
	return comovingElement(cosmo, z.value)
end
function comovingElement(cosmo::CosmologicalModel, a::ScaleFactor)
	return comovingElement(cosmo, convert(Redshift, a))
end


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