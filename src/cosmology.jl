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
	function CosmologicalModel{C}(cosmo::C; z = nothing) where {C}
		if isnothing(z)
			z1 = ScaleLogarithmicNegativeRange(-1.0 + 1e-10, -0.01, 91)
			z2 = ScaleLinearRange(-0.01, 0.01, 21)
			z3 = ScaleLogarithmicRange(0.01, 1e3, 51)
			z = vcat(z1, z2, z3)
			unique!(z)
		end
	
		# distance arrays
		dC, dL, dP = [], [], []

		@simd for i in eachindex(z)
			dC0 = comoving_radial_dist(cosmo, 0., z[i]) |> u"m"
			dL0 = luminosity_dist(cosmo, z[i]) |> u"m"
			dP0 = (lookback_time(cosmo, z[i]) |> u"s") * SpeedOfLightInVacuum |> u"m"
			@inbounds push!(dC, dC0)
			@inbounds push!(dL, dL0)
			@inbounds push!(dP, dP0)
		end
	
		# get indices that sort arrays, as required by Interpolations.jl
		idxC = sortperm(dC)
		idxL = sortperm(dL)
		idxP = sortperm(dP)
	
		# functions 
		c2z = LinearInterpolation(ustrip.(dC[idxC]), z[idxC])
		l2z = LinearInterpolation(ustrip.(dL[idxL]), z[idxL])
		p2z = LinearInterpolation(ustrip.(dP[idxP]), z[idxP])


		return new{C}(cosmo, c2z, l2z, p2z)
	end
end

# convenience constructors
const CosmologyPlanck(; z = nothing) = CosmologicalModel{Cosmology.FlatLCDM{Float64}}(cosmology(); z = z)



# ---------------------------------------------------------------------------------- #
# 
# @doc """
# 	comovingDistanceToRedshift(cosmo, d)
# 	comovingDistanceToRedshift(cosmo, d1, d2)

# Convert comoving distance to redshift
# """
redshiftsToComovingDistance(cosmo::CosmologicalModel, z0::Float64, z1::Float64) = comoving_radial_dist(cosmo.cosmology, z0, z1) 
redshiftToComovingDistance(cosmo::CosmologicalModel, z0::Float64) = redshiftsToComovingDistance(cosmo, 0., z0)

redshiftToLuminosityDistance(cosmo::CosmologicalModel, z0::Float64) = luminosity_dist(cosmo.cosmology, z0)
redshifts2ToLuminosityDistance(cosmo::CosmologicalModel, z0::Float64, z1::Float64) = redshiftToLuminosityDistance(cosmo, z0) - redshiftToLuminosityDistance(cosmo, z1)

redshiftToLookbackTime(cosmo::CosmologicalModel, z0::Float64) = lookback_time(cosmo.cosmology, z0)

redshiftToLightTravelDistance(cosmo::CosmologicalModel, z0::Float64) = ((redshiftToLookbackTime(cosmo, z0) |> u"s") * SpeedOfLightInVacuum) |> u"Mpc"
redshiftsToLightTravelDistance(cosmo::CosmologicalModel, z0::Float64, z1::Float64) = redshiftToLightTravelDistance(cosmo, z0) - redshiftToLightTravelDistance(cosmo, z1)

comovingDistanceToRedshift(cosmo::CosmologicalModel, d0::Float64) = cosmo._comovingDistance2Redshift(d0)
comovingDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.AbstractQuantity) = (dimension(d0) == 𝐋) ? comovingDistanceToRedshift(cosmo, ustrip(d0 |> u"m")) : throw(DimensionMismatch("Dimension of provided quantity is not distance."))

lightTravelDistanceToRedshift(cosmo::CosmologicalModel, d0::Float64) = cosmo._lightTravelDistance2Redshift(d0)
lightTravelDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.AbstractQuantity) = (dimension(d0) == 𝐋) ? lightTravelDistanceToRedshift(cosmo, ustrip(d0 |> u"m")) : throw(DimensionMismatch("Dimension of provided quantity is not distance."))

luminosityDistanceToRedshift(cosmo::CosmologicalModel, d0::Float64) = cosmo._luminosityDistance2Redshift(d0)
luminosityDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.AbstractQuantity) = (dimension(d0) == 𝐋) ? luminosityDistanceToRedshift(cosmo, ustrip(d0 |> u"m")) : throw(DimensionMismatch("Dimension of provided quantity is not distance."))


# ---------------------------------------------------------------------------------- #
# 