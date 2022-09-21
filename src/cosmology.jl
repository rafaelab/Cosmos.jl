# ---------------------------------------------------------------------------------- #
#
@doc """
General struct to hold a `Cosmology.jl` object of type `AbstractCosmology`.
The only member of this struct that should be accessed is cosmology. 
The others are only useful for the functions.

Any cosmology can be built using the following parameters:
. `h`: dimensionless Hubble constant \\
. `Ωm`: matter density \\
. `Ωr`: radiation density \\
. `Ωk`: curvature density \\
. `ΩΛ`: dark energy density \\
. `wEOSΛ`: tuple with parameters of the equation of state for dark energy: `w = w_0 + w_a (1 - a)` \\

The default constructors can be built using only the first 3 or 4 parameters.
"""
struct CosmologicalModel{C <: Cosmology.AbstractCosmology, T <: Real}
	h::T
	ΩΛ::T
	Ωm::T
	Ωr::T
	Ωk::T
	wEOSΛ::Tuple{T, T}
	_cosmology::C
	_comovingDistance2Redshift
	_luminosityDistance2Redshift
	_lightTravelDistance2Redshift
	_comovingTransverse2Redshift
	_angularDiameterDistance2Redshift
	function CosmologicalModel{C, T}(h::Real, Ωm::Real, Ωr::Real, Ωk::Real, wEOSΛ::Tuple; z::Union{Nothing, Vector{Z}} = nothing) where {C <: Cosmology.AbstractCosmology, T <: Real, Z <: Real}
		cosmo = cosmology(; h = h, OmegaK = Ωk, OmegaM = Ωm, OmegaR = Ωr, w0 = wEOSΛ[1], wa = wEOSΛ[2])
		ΩΛ = cosmo.Ω_Λ

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
		dC, dL, dP, dT, dA = T[], T[], T[], T[], T[]

		@simd for i in eachindex(z)
			dC0 = comoving_radial_dist(cosmo, 0., z[i]) |> u"m"
			dL0 = luminosity_dist(cosmo, z[i]) |> u"m"
			dP0 = (lookback_time(cosmo, z[i]) |> u"s") * SpeedOfLightInVacuum |> u"m"
			dT0 = comoving_transverse_dist(cosmo, z[i]) |> u"m"
			dA0 = angular_diameter_dist(cosmo, z[i]) |> u"m"
			@inbounds push!(dC, ustrip(dC0 |> u"m"))
			@inbounds push!(dL, ustrip(dL0 |> u"m"))
			@inbounds push!(dP, ustrip(dP0 |> u"m"))
			@inbounds push!(dT, ustrip(dT0 |> u"m"))
			@inbounds push!(dA, ustrip(dA0 |> u"m"))
		end
	
		# get indices that sort arrays, as required by Interpolations.jl
		idxC = sortperm(dC)
		idxL = sortperm(dL)
		idxP = sortperm(dP)
		idxT = sortperm(dT)
		idxA = sortperm(dA)
	
		# functions: for `Interpolations` version < 0.14
		# c2z = LinearInterpolation(ustrip.(dC[idxC]), z[idxC])
		# l2z = LinearInterpolation(ustrip.(dL[idxL]), z[idxL])
		# p2z = LinearInterpolation(ustrip.(dP[idxP]), z[idxP])

		c2z = linear_interpolation(dC[idxC], z[idxC])
		l2z = linear_interpolation(dL[idxL], z[idxL])
		p2z = linear_interpolation(dP[idxP], z[idxP])
		t2z = linear_interpolation(dT[idxT], z[idxP])
		a2z = linear_interpolation(dA[idxA], z[idxP])

		return new{typeof(cosmo), typeof(h)}(h, ΩΛ, Ωm, Ωr, Ωk, wEOSΛ, cosmo, c2z, l2z, p2z, t2z, a2z)
	end
end

function CosmologicalModel(h::Real, Ωm::Real, Ωr::Real, Ωk::Real, wEOSΛ::Tuple; z::Union{Nothing, Vector{Z}} = nothing) where {Z <: Real}
	T = promote_type(typeof(h), typeof(Ωm), typeof(Ωk), typeof(Ωr), typeof(wEOSΛ[1]), typeof(wEOSΛ[2]))
	h = T(h)
	Ωk = T(Ωk)
	Ωr = T(Ωr) 
	wEOSΛ = (T(first(wEOSΛ)), T(last(wEOSΛ)))

	cosmo = cosmology(; h = h, OmegaK = Ωk, OmegaM = Ωm, OmegaR = Ωr, w0 = wEOSΛ[1], wa = wEOSΛ[2])
	C = typeof(cosmo)

	return CosmologicalModel{C, T}(h, Ωm, Ωr, Ωk, wEOSΛ; z = z)
end

CosmologicalModel(h::Real, Ωm::Real, Ωr::Real, Ωk::Real; z::Union{Nothing, Vector{Z}} = nothing) where {Z <: Real} = CosmologicalModel(h, Ωm, Ωr, Ωk, (-1., 0.); z = z) 

CosmologicalModel(h::Real, Ωm::Real, Ωr::Real; z::Union{Nothing, Vector{Z}} = nothing) where {Z <: Real} = CosmologicalModel(h, Ωm, Ωr, 0., (-1., 0.); z = z) 

CosmologicalModel(cosmo::FlatLCDM{T}; z::Union{Nothing, Vector{Z}} = nothing) where {T <: Real, Z <: Real} = CosmologicalModel(cosmo.h, cosmo.Ω_m, cosmo.Ω_r, 0., (-1., 0.); z = z)

CosmologicalModel(cosmo::Union{OpenLCDM{T}, ClosedLCDM{T}}; z::Union{Nothing, Vector{Z}} = nothing) where {T <: Real, Z <: Real} = CosmologicalModel(cosmo.h, cosmo.Ω_m, cosmo.Ω_r, cosmo.Ω_k, (-1., 0.); z = z)

CosmologicalModel(cosmo::FlatWCDM{T}; z::Union{Nothing, Vector{Z}} = nothing) where {T <: Real, Z <: Real} = CosmologicalModel(cosmo.h, cosmo.Ω_m, cosmo.Ω_r, 0., (cosmo.w0, cosmo.wa); z = z)

CosmologicalModel(cosmo::Union{ClosedWCDM{T}, OpenWCDM{T}}; z::Union{Nothing, Vector{Z}} = nothing) where {T <: Real, Z <: Real} = CosmologicalModel(cosmo.h, cosmo.Ω_m, cosmo.Ω_r, cosmo.Ω_k, (cosmo.w0, cosmo.wa); z = z)



		
# ---------------------------------------------------------------------------------- #
# 
@doc """
Latest Planck's cosmology.
"""
const CosmologyPlanck(; z::Union{Nothing, Vector{Z}} = nothing) where {Z <: Real} = CosmologicalModel(cosmology(); z = z)


# ---------------------------------------------------------------------------------- #
# 
@doc """ 
	dimensionlessHubbleParameter(cosmology, z)

Compute the dimensionless Hubble parameter (often dubbed `E`):
  ``E(z) = \\frac{H(z)}{H(0)} \\sqrt{\\Omega_r (1 + z)^4 + \\Omega_m (1 + z)^3 + \\Omega_k (1 + z)^2 + \\Omega_\\Lambda}``
This follows the definition from Peebles 1993 (p. ~310-322), adopted by Hogg, arXiv:astro-ph/9905116.
"""
dimensionlessHubbleParameter(cosmo::CosmologicalModel, z::Real) = E(cosmo._cosmology, z)


# ---------------------------------------------------------------------------------- #
# 
@doc """ 
	hubbleParameter(cosmology, z)

Compute the Hubble parameter H(z).
"""
hubbleParameter(cosmo::CosmologicalModel, z::Real) = H(cosmo._cosmology, z)


# ---------------------------------------------------------------------------------- #
# 
@doc """ 
	scaleFactor(z)

Compute the scale factor for a given cosmology. 
This should be the same for all cosmologies, by definition, but passing this argument fixes the correct type.
"""
scaleFactor(cosmo::CosmologicalModel, z::Real) = eltype(cosmo)(1 / (1 + z))
scaleFactor(z::Real) = promote_type(Float64, typeof(z))(1. / (1. + z))


# ---------------------------------------------------------------------------------- #
# 
@doc """ 
	hubbleDistance(cosmo)
	hubbleDistance(cosmo, z)

Computes the Hubble distance for a given cosmology and possibly at a given redshift.
"""
hubbleDistance(cosmo::CosmologicalModel) = hubble_dist0(cosmo._cosmology)
hubbleDistance(cosmo::CosmologicalModel, z::Real) = hubble_dist(cosmo.cosmology, z)


# ---------------------------------------------------------------------------------- #
# 
@doc """ 
	hubbleTime(cosmo)
	hubbleTime(cosmo, z)

Computes the Hubble distance for a given cosmology and possibly at a given redshift.
"""
hubbleTime(cosmo::CosmologicalModel) = hubble_time0(cosmo._cosmology)
hubbleTime(cosmo::CosmologicalModel, z::Real) = hubble_time(cosmo._cosmology, z)


# ---------------------------------------------------------------------------------- #
# 
@doc """ 
	ageOfUniverse(cosmo)
	ageOfUniverse(cosmo, z)

Computes the Hubble distance for a given cosmology and possibly at a given redshift.
"""
ageOfUniverse(cosmo::CosmologicalModel) = age(cosmo._cosmology)
ageOfUniverse(cosmo::CosmologicalModel, z::Real) = age(cosmo._cosmology, z)


# ---------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToComovingDistance(cosmology, z)
	redshiftToComovingDistance(cosmology, z1, z2)

Computes the comoving distance between two redshifts.
If a second redshift value is not provided, it is assumed to be 0.
"""
redshiftToComovingDistance(cosmo::CosmologicalModel, z0::Real, z1::Real) = comoving_radial_dist(cosmo._cosmology, z0, z1) 
redshiftToComovingDistance(cosmo::CosmologicalModel, z::Real) = comoving_radial_dist(cosmo._cosmology, zero(eltype(cosmo)), z)


# ---------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToLuminosityDistance(cosmology, z)
	redshiftToLuminosityDistance(cosmology, z1, z2)

Computes the luminosity distance between two redshifts.
If a second redshift value is not provided, it is assumed to be 0.
"""
redshiftToLuminosityDistance(cosmo::CosmologicalModel, z0::Real) = luminosity_dist(cosmo._cosmology, z0)
redshiftToLuminosityDistance(cosmo::CosmologicalModel, z0::Real, z1::Real) = redshiftToLuminosityDistance(cosmo, z0) - redshiftToLuminosityDistance(cosmo, z1)


# ---------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToLookbackTime(cosmology, z)

Computes the lookback_time corresponding to a given redshift.
"""
redshiftToLookbackTime(cosmo::CosmologicalModel, z0::Real) = lookback_time(cosmo._cosmology, z0)


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
	redshiftToAngularDiameterDistance(cosmology, z)

Computes the angular diameter distance at a given redshift.
The angular diameter distance takes into account the object's transverse and its angular size.
"""
redshiftToAngularDiameterDistance(cosmo::CosmologicalModel, z0::Real) = angular_diameter_dist(cosmo._cosmology, z0) |> u"Mpc"
redshiftToAngularDiameterDistance(cosmo::CosmologicalModel, z1::Real, z2::Real) = angular_diameter_dist(cosmo._cosmology, z1, z2) |> u"Mpc"


# ---------------------------------------------------------------------------------- #
# 
@doc """
	redshiftToTransverseComovingDistance(cosmology, z)

Computes the angular diameter distance at a given redshift.
The angular diameter distance takes into account the object's transverse and its angular size.
"""
redshiftToTransverseComovingDistance(cosmo::CosmologicalModel, z1::Real, z2::Real) = transverse_comoving_dist(cosmo._cosmology, z1, z2) |> u"Mpc"
redshiftToTransverseComovingDistance(cosmo::CosmologicalModel, z::Real) = redshiftToTransverseComovingDistance(cosmo, zero(eltype(cosmo)), z)


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
@doc """
	comovingTransverseDistanceToRedshift(cosmology, z)
	comovingTransverseDistanceToRedshift(cosmology, z1, z2)

Computes the redshift corresponding to a given transverse comoving distance.
"""
comovingTransverseDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._comovingTransverseDistance2Redshift(d0)
comovingTransverseDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.AbstractQuantity) = isLengthDimension(d0) && comovingTransverseDistanceToRedshift(cosmo, ustrip(d0 |> u"m"))


# ---------------------------------------------------------------------------------- #
# 
@doc """
	angularDiameterDistanceToRedshift(cosmology, z)
	angularDiameterDistanceToRedshift(cosmology, z1, z2)

Computes the redshift corresponding to a given angular diameter distance.
"""
angularDiameterDistanceToRedshift(cosmo::CosmologicalModel, d0::Real) = cosmo._angularDiameterDistance2Redshift(d0)
angularDiameterDistanceToRedshift(cosmo::CosmologicalModel, d0::Unitful.AbstractQuantity) = isLengthDimension(d0) && angularDiameterDistanceToRedshift(cosmo, ustrip(d0 |> u"m"))


# ---------------------------------------------------------------------------------- #
# 
@doc """
	comovingVolume(cosmology, z)

Calculates the comoving volume at a given redshift.
"""
comovingVolume(cosmo::CosmologicalModel, z::Real) = comoving_volume(cosmo._cosmology, z)


# ---------------------------------------------------------------------------------- #
# 
@doc """
	comovingVolumeElement(cosmology, z)

Calculates the comoving volume element at a given redshift.
"""
comovingVolumeElement(cosmo::CosmologicalModel, d0::Real) = comoving_volume_element(cosmo._cosmology, z)

# ---------------------------------------------------------------------------------- #
# 
@doc """
	calculateDensityNeutrinos(cosmology, Nν, Tcmb)

Compute the density of relativistic components.

# Input
. `Nν`: effective number of neutrino species (3.04 according to Planck) \\
. `Tcmb`: current CMB temperature \\
"""
calculateDensityNeutrinos(cosmo::CosmologicalModel, Nν::Real, Tcmb::Real) = eltype(cosmo)((7. / 8.) * (4. / 11.) ^ (4 / 3) * Nν * calculateDensityPhotons(cosmo, Tcmb))
	


# ---------------------------------------------------------------------------------- #
# 
@doc """
	calculateDensityPhotons(cosmology, Nν, Tcmb)

Compute the density of photons.

# Input
. `Tcmb`: current CMB temperature \\
"""
calculateDensityPhotons(cosmo::CosmologicalModel, Tcmb::Real) = eltype(cosmo)(4.4813e-7 * (Tcmb ^ 2 / cosmo.h) ^ 2)


# ---------------------------------------------------------------------------------- #
# 
@doc """
	eltype(cosmology)

Overload Base functions. 
The returned type is essentially the data type of the underlying cosmology from `Cosmology.jl` (Float64, by default).
"""
Base.eltype(cosmo::CosmologicalModel) = typeof(cosmo._cosmology.h)


# ---------------------------------------------------------------------------------- #
# 
# check if dimension provided is correct
@inline function isLengthDimension(d::Unitful.AbstractQuantity) 
	if dimension(d0) ≠ 𝐋
		throw(DimensionMismatch("Dimension of provided quantity is not distance."))
	end

	return true
end

# ---------------------------------------------------------------------------------- #
# 