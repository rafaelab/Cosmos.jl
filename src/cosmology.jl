# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	dimensionlessHubbleParameter(cosmology, z)

Compute the dimensionless Hubble parameter (often dubbed `E`):
  ``E(z) = \\frac{H(z)}{H(0)} \\sqrt{\\Omega_r (1 + z)^4 + \\Omega_m (1 + z)^3 + \\Omega_k (1 + z)^2 + \\Omega_\\Lambda}``
This follows the definition from Peebles 1993 (p. ~310-322), adopted by Hogg, arXiv:astro-ph/9905116.
"""
dimensionlessHubbleParameter(cosmo::CosmologicalModel, z::Real) = E(cosmo._cosmology, z)


# ----------------------------------------------------------------------------------------------- #
@doc """ 
	hubbleParameter(cosmology, z)

Compute the Hubble parameter H(z).
"""
hubbleParameter(cosmo::CosmologicalModel, z::Real) = H(cosmo._cosmology, z)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	hubbleConstant(cosmology, z)

Compute the present-day Hubble constant.
"""
hubbleConstant(cosmo::CosmologicalModel) = hubbleParameter(cosmo, 0.)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	scaleFactor(z)

Compute the scale factor for a given cosmology. 
This should be the same for all cosmologies, by definition, but passing this argument fixes the correct type.
"""
scaleFactor(cosmo::CosmologicalModel, z::Real) = eltype(cosmo)(1 / (1 + z))
scaleFactor(z::Real) = promote_type(Float64, typeof(z))(1. / (1. + z))


# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	hubbleDistance(cosmo)
	hubbleDistance(cosmo, z)

Computes the Hubble distance for a given cosmology and possibly at a given redshift.
"""
hubbleDistance(cosmo::CosmologicalModel) = hubble_dist0(cosmo._cosmology)
hubbleDistance(cosmo::CosmologicalModel, z::Real) = hubble_dist(cosmo._cosmology, z)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	hubbleTime(cosmo)
	hubbleTime(cosmo, z)

Computes the Hubble distance for a given cosmology and possibly at a given redshift.
"""
hubbleTime(cosmo::CosmologicalModel) = hubble_time0(cosmo._cosmology)
hubbleTime(cosmo::CosmologicalModel, z::Real) = hubble_time(cosmo._cosmology, z)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """ 
	ageOfUniverse(cosmo)
	ageOfUniverse(cosmo, z)

Computes the Hubble distance for a given cosmology and possibly at a given redshift.
"""
ageOfUniverse(cosmo::CosmologicalModel) = age(cosmo._cosmology)
ageOfUniverse(cosmo::CosmologicalModel, z::Real) = age(cosmo._cosmology, z)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	comovingVolume(cosmology, z)

Calculates the comoving volume at a given redshift.
"""
comovingVolume(cosmo::CosmologicalModel, z::Real) = comoving_volume(cosmo._cosmology, z)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	comovingVolumeElement(cosmology, z)

Calculates the comoving volume element at a given redshift.
"""
comovingVolumeElement(cosmo::CosmologicalModel, z0::Real) = comoving_volume_element(cosmo._cosmology, z0)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	comovingElement(cosmology, z)

Calculates the comoving line element at a given redshift.
This is Hogg's eq. 28 adjusted.
NOTE: check nomenclature
"""
comovingElement(cosmo::CosmologicalModel, z0::Float64) = hubbleDistance(cosmo, z0)  / E(cosmo._cosmology, z0) / (1 + z0)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	calculateDensityNeutrinos(cosmology, Nν, Tcmb)

Compute the density of relativistic components.

# Input
. `Nν`: effective number of neutrino species (3.04 according to Planck) \\
. `Tcmb`: current CMB temperature \\
"""
calculateDensityNeutrinos(cosmo::CosmologicalModel, Nν::Real, Tcmb::Real) = eltype(cosmo)((7. / 8.) * (4. / 11.) ^ (4 / 3) * Nν * calculateDensityPhotons(cosmo, Tcmb))
	

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	calculateDensityPhotons(cosmology, Nν, Tcmb)

Compute the density of photons.

# Input
. `Tcmb`: current CMB temperature \\
"""
calculateDensityPhotons(cosmo::CosmologicalModel, Tcmb::Real) = eltype(cosmo)(4.4813e-7 * (Tcmb ^ 2 / cosmo.h) ^ 2)


# ----------------------------------------------------------------------------------------------- #
# 