# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	computeCriticalDensity(cosmology; z = 0.)

Computes the critical density of the universe at a given redshift:
  ρc = 3 H^2(z) / 8π G

# Input 
. `cosmology::CosmologicalModel`: the cosmological model of interest \\
. `z::Real`: the redshift at which to compute the density \\
"""
function computeCriticalDensity(cosmology::CosmologicalModel; z::Real = 0.)
	return upreferred(3. * hubbleParameter(cosmology, z) ^ 2 / 8π / G)
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	computeMatterDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
  ρm = Ωm * ρc

# Input 
. `cosmology::CosmologicalModel`: the cosmological model of interest \\
. `z::Real`: the redshift at which to compute the density \\
"""
function computeMatterDensity(cosmology::CosmologicalModel; z::Real = 0.)
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.Ωm * ρc)
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	computeRadiationDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
  ρr = Ωr * ρc

# Input 
. `cosmology::CosmologicalModel`: the cosmological model of interest \\
. `z::Real`: the redshift at which to compute the density \\
"""
function computeRadiationDensity(cosmology::CosmologicalModel; z::Real = 0.)
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.Ωr * ρc)
end

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	computeCurvatureDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
  ρk = Ωk * ρc

# Input 
. `cosmology::CosmologicalModel`: the cosmological model of interest \\
. `z::Real`: the redshift at which to compute the density \\
"""
function computeCurvatureDensity(cosmology::CosmologicalModel; z::Real = 0.)
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.Ωk * ρc)
end

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	computeDarkEnergyDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
  ρΛ = ΩΛ * ρc

# Input 
. `cosmology::CosmologicalModel`: the cosmological model of interest \\
. `z::Real`: the redshift at which to compute the density \\
"""
function computeDarkEnergyDensity(cosmology::CosmologicalModel; z::Real = 0.)
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.ΩΛ * ρc)
end

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	computeBaryonDensity(cosmology; z = 0.)

Computes the density of matter in the universe at a given redshift:
  ρb = Ωb * ρc

# Input 
. `cosmology::CosmologicalModel`: the cosmological model of interest \\
. `z::Real`: the redshift at which to compute the density \\
"""
function computeBaryonDensity(cosmology::CosmologicalModel; z::Real = 0.)
	if cosmology.Ωb < 0
		throw(ArgumentError("Cannot compute the baryon density because the baryon fraction was not provided to the `CosmologicalModel`."))
	end
	ρc = computeCriticalDensity(cosmology; z = z)
	return upreferred(cosmology.Ωb * ρc)
end


# ----------------------------------------------------------------------------------------------- #
# 