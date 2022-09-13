# ---------------------------------------------------------------------------------- #
#
@doc """
General struct to hold a `Cosmology.jl` object of type `AbstractCosmology`.
The only member of this struct that should be accessed is cosmology. 
The others are only useful for the functions.
"""
struct CosmologicalModel{C <: Cosmology.AbstractCosmology}
	cosmology::C
	_comovingDistance2Redshift
	_luminosityDistance2Redshift
	_lightTravelDistance2Redshift
	function CosmologicalModel{C}(cosmo::C; z::Union{Nothing, Vector{Z}} = nothing) where {C <: Cosmology.AbstractCosmology, Z <: Real}
		# infer general type
		T = typeof(cosmo.h)

		# prepare redshifts
		if isnothing(z)
			z = T[]
			append!(z, -10. .^ collect(range(-3., 0.; length = 91)))
			append!(z, collect(range(-0.001, 0.001; length = 31)))
			append!(z, 10 .^ collect(range(-3., 3.; length = 61)))
			unique!(z)
		else
			z = convert(Vector{T}, z)
		end
		z = z[z .> -1.0]
		sort!(z)

	
		# distance arrays
		dC, dL, dP = T[], T[], T[]

		@simd for i in eachindex(z)
			dC0 = comoving_radial_dist(cosmo, 0., z[i]) |> u"m"
			dL0 = luminosity_dist(cosmo, z[i]) |> u"m"
			dP0 = (lookback_time(cosmo, z[i]) |> u"s") * SpeedOfLightInVacuum |> u"m"
			@inbounds push!(dC, ustrip(dC0 |> u"m"))
			@inbounds push!(dL, ustrip(dL0 |> u"m"))
			@inbounds push!(dP, ustrip(dP0 |> u"m"))
		end
	
		# get indices that sort arrays, as required by Interpolations.jl
		idxC = sortperm(dC)
		idxL = sortperm(dL)
		idxP = sortperm(dP)
	
		# functions: for `Interpolations` version < 0.14
		# c2z = LinearInterpolation(ustrip.(dC[idxC]), z[idxC])
		# l2z = LinearInterpolation(ustrip.(dL[idxL]), z[idxL])
		# p2z = LinearInterpolation(ustrip.(dP[idxP]), z[idxP])

		println("dC, ", dC[idxC])

		c2z = linear_interpolation(dC[idxC], z[idxC])
		l2z = linear_interpolation(dL[idxL], z[idxL])
		p2z = linear_interpolation(dP[idxP], z[idxP])


		return new{C}(cosmo, c2z, l2z, p2z)
	end
end

# ---------------------------------------------------------------------------------- #
# 
@doc """
Define Planck's cosmology.
"""
const CosmologyPlanck(; z = nothing) = CosmologicalModel{Cosmology.FlatLCDM{Float64}}(cosmology(); z = z)

# ---------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToComovingDistance(cosmology, z)
	redshiftToComovingDistance(cosmology, z1, z2)

Computes the comoving distance between two redshifts.
If a second redshift value is not provided, it is assumed to be 0.
"""
redshiftToComovingDistance(cosmo::CosmologicalModel, z0::Real, z1::Real) = comoving_radial_dist(cosmo.cosmology, z0, z1) 
redshiftToComovingDistance(cosmo::CosmologicalModel, z::Real) = comoving_radial_dist(cosmo.cosmology, zero(eltype(cosmo)), z)

# ---------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToLuminosityDistance(cosmology, z)
	redshiftToLuminosityDistance(cosmology, z1, z2)

Computes the luminosity distance between two redshifts.
If a second redshift value is not provided, it is assumed to be 0.
"""
redshiftToLuminosityDistance(cosmo::CosmologicalModel, z0::Real) = luminosity_dist(cosmo.cosmology, z0)
redshiftToLuminosityDistance(cosmo::CosmologicalModel, z0::Real, z1::Real) = redshiftToLuminosityDistance(cosmo, z0) - redshiftToLuminosityDistance(cosmo, z1)


# ---------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToLookbackTime(cosmology, z)

Computes the lookback_time corresponding to a given redshift.
"""
redshiftToLookbackTime(cosmo::CosmologicalModel, z0::Real) = lookback_time(cosmo.cosmology, z0)


# ---------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToLightTravelDistance(cosmology, z)

Computes the light-travel distance between two redshifts.
"""
redshiftToLightTravelDistance(cosmo::CosmologicalModel, z0::Real) = (redshiftToLookbackTime(cosmo, z0) * SpeedOfLightInVacuum) |> u"Mpc"


# ---------------------------------------------------------------------------------- #
# 
@doc """
	comovingDistanceToRedshift(cosmology, z)
	comovingDistanceToRedshift(cosmology, z1, z2)

Computes the redshift corresponding to a given comoving distance.
"""
comovingDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._comovingDistance2Redshift(d0)
comovingDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.AbstractQuantity) = isLengthDimension(d0) && comovingDistanceToRedshift(cosmo, ustrip(d0 |> u"m"))


# ---------------------------------------------------------------------------------- #
# 
@doc """
	lightTravelDistanceToRedshift(cosmology, z)

Computes the redshift corresponding to a given light-travel distance.
"""
lightTravelDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._lightTravelDistance2Redshift(d0)
lightTravelDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.AbstractQuantity) = isLengthDimension(d0) && lightTravelDistanceToRedshift(cosmo, ustrip(d0 |> u"m")) 


# ---------------------------------------------------------------------------------- #
# 
@doc """
	luminosityDistanceToRedshift(cosmology, z)

Computes the redshift corresponding to a given luminosity distance.
"""
luminosityDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._luminosityDistance2Redshift(d0)
luminosityDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.AbstractQuantity) = isLengthDimension(d0) &&  luminosityDistanceToRedshift(cosmo, ustrip(d0 |> u"m")) 


# ---------------------------------------------------------------------------------- #
# 
# Overload Base functions
Base.eltype(cosmo::CosmologicalModel) = typeof(cosmo.cosmology.h)

# ---------------------------------------------------------------------------------- #
# 
# check if dimension provided is correct
@inline function isLengthDimension(d::Unitful.AbstractQuantity) 
	if dimension(d0) ≠ 𝐋
		throw(DimensionMismatch("Dimension of provided quantity is not distance."))
	end

	return true
end
_

# ---------------------------------------------------------------------------------- #
# 