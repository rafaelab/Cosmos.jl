# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToComovingDistance(cosmology, z)
	redshiftToComovingDistance(cosmology, z1, z2)

Computes the comoving distance between two redshifts.
If a second redshift value is not provided, it is assumed to be 0.
"""
redshiftToComovingDistance(cosmo::CosmologicalModel, z0::Real, z1::Real) = comoving_radial_dist(cosmo._cosmology, z0, z1) 
redshiftToComovingDistance(cosmo::CosmologicalModel, z::Real) = comoving_radial_dist(cosmo._cosmology, zero(eltype(cosmo)), z)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToLuminosityDistance(cosmology, z)
	redshiftToLuminosityDistance(cosmology, z1, z2)

Computes the luminosity distance between two redshifts.
If a second redshift value is not provided, it is assumed to be 0.
"""
redshiftToLuminosityDistance(cosmo::CosmologicalModel, z0::Real) = luminosity_dist(cosmo._cosmology, z0)
redshiftToLuminosityDistance(cosmo::CosmologicalModel, z0::Real, z1::Real) = redshiftToLuminosityDistance(cosmo, z0) - redshiftToLuminosityDistance(cosmo, z1)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToLookbackTime(cosmology, z)

Computes the lookback_time corresponding to a given redshift.
"""
redshiftToLookbackTime(cosmo::CosmologicalModel, z0::Real) = lookback_time(cosmo._cosmology, z0)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToLightTravelDistance(cosmology, z)

Computes the light-travel distance between two redshifts.
"""
redshiftToLightTravelDistance(cosmo::CosmologicalModel, z0::Real) = (redshiftToLookbackTime(cosmo, z0) * SpeedOfLightInVacuum) |> u"Mpc"


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToAngularDiameterDistance(cosmology, z)

Computes the angular diameter distance at a given redshift.
The angular diameter distance takes into account the object's transverse and its angular size.
"""
redshiftToAngularDiameterDistance(cosmo::CosmologicalModel, z0::Real) = angular_diameter_dist(cosmo._cosmology, z0) |> u"Mpc"
redshiftToAngularDiameterDistance(cosmo::CosmologicalModel, z1::Real, z2::Real) = angular_diameter_dist(cosmo._cosmology, z1, z2) |> u"Mpc"


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToTransverseComovingDistance(cosmology, z)

Computes the angular diameter distance at a given redshift.
The angular diameter distance takes into account the object's transverse and its angular size.
"""
redshiftToTransverseComovingDistance(cosmo::CosmologicalModel, z1::Real, z2::Real) = transverse_comoving_dist(cosmo._cosmology, z1, z2) |> u"Mpc"
redshiftToTransverseComovingDistance(cosmo::CosmologicalModel, z::Real) = redshiftToTransverseComovingDistance(cosmo, zero(eltype(cosmo)), z)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	comovingDistanceToRedshift(cosmology, z)
	comovingDistanceToRedshift(cosmology, z1, z2)

Computes the redshift corresponding to a given comoving distance.
"""
comovingDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._comovingDistance2Redshift(d0)
comovingDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.Length) = comovingDistanceToRedshift(cosmo, ustrip(d0 |> u"m"))


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	lightTravelDistanceToRedshift(cosmology, z)

Computes the redshift corresponding to a given light-travel distance.
"""
lightTravelDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._lightTravelDistance2Redshift(d0)
lightTravelDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.Length) = lightTravelDistanceToRedshift(cosmo, ustrip(d0 |> u"m")) 


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	luminosityDistanceToRedshift(cosmology, z)

Computes the redshift corresponding to a given luminosity distance.
"""
luminosityDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._luminosityDistance2Redshift(d0)
luminosityDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.Length) = luminosityDistanceToRedshift(cosmo, ustrip(d0 |> u"m")) 


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	comovingTransverseDistanceToRedshift(cosmology, z)

Computes the redshift corresponding to a given transverse comoving distance.
"""
comovingTransverseDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._comovingTransverseDistance2Redshift(d0)
comovingTransverseDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.Length) = comovingTransverseDistanceToRedshift(cosmo, ustrip(d0 |> u"m"))


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	angularDiameterDistanceToRedshift(cosmology, z)

Computes the redshift corresponding to a given angular diameter distance.
"""
angularDiameterDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._angularDiameterDistance2Redshift(d0)
angularDiameterDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.Length) = angularDiameterDistanceToRedshift(cosmo, ustrip(d0 |> u"m"))


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	lightTravelToComovingDistance(cosmology, d)

Computes the comoving distance corresponding to a given light-travel distance.
"""
lightTravelToComovingDistance(cosmo::CosmologicalModel, d0::Real) = redshiftToComovingDistance(cosmo, cosmo._lightTravelDistance2Redshift(d0))
lightTravelToComovingDistance(cosmo::CosmologicalModel, d0::Unitful.Length) = lightTravelToComovingDistance(cosmo, ustrip(d0 |> u"m"))

# ----------------------------------------------------------------------------------------------- #