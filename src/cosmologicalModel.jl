# ----------------------------------------------------------------------------------------------- #
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


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Latest Planck's cosmology.
"""
const CosmologyPlanck(; z::Union{Nothing, Vector{Z}} = nothing) where {Z <: Real} = CosmologicalModel(cosmology(); z = z)


# ----------------------------------------------------------------------------------------------- #
# 
# Overload Base functions. 
# The returned type is essentially the data type of the underlying cosmology from `Cosmology.jl` (Float64, by default).
Base.eltype(cosmo::CosmologicalModel) = typeof(cosmo._cosmology.h)



# ----------------------------------------------------------------------------------------------- #