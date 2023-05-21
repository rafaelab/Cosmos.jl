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
	function CosmologicalModel{C, T}(cosmo::C; z::Union{Nothing, Vector{Z}} = nothing) where {C <: Cosmology.AbstractCosmology, T <: Real, Z <: Real}		
		h = cosmo.h
		Ωk = (C isa Union{Cosmology.AbstractClosedCosmology, Cosmology.AbstractOpenCosmology}) ? cosmo.Ω_k : zero(T)
		Ωr = cosmo.Ω_r
		Ωm = cosmo.Ω_m
		ΩΛ = cosmo.Ω_Λ
		wEOSΛ = (C isa Union{FlatWCDM, OpenWCDM, ClosedWCDM}) ? (cosmo.w0, cosmo.wa) : (- one(T), zero(T))
		
		# prepare redshifts
		if isnothing(z)
			z = T[]
			append!(z, -10. .^ collect(range(-3., 0.; length = 91)))
			append!(z, collect(range(-0.001, 0.001; length = 31)))
			append!(z, collect(range(0.001, 10.; length = 81)))
			append!(z, collect(range(1., 2.; length = 21)))
			append!(z, 10 .^ collect(range(2., 4.; length = 41)))
			unique!(z)
		else
			z = convert(Vector{T}, z)
		end
		z = z[z .> -1.0]
		sort!(z)

		# distance arrays
		dC, dL, dP, dT, dA = T[], T[], T[], T[], T[]

		@simd for i in eachindex(z)
			dC0 = comoving_radial_dist(cosmo, z[i]) |> u"m"
			dL0 = luminosity_dist(cosmo, z[i]) |> u"m"
			dP0 = (lookback_time(cosmo, z[i]) |> u"s") * SpeedOfLightInVacuum |> u"m"
			dT0 = comoving_transverse_dist(cosmo, z[i]) |> u"m"
			dA0 = angular_diameter_dist(cosmo, z[i]) |> u"m"
			@inbounds push!(dC, ustrip.(dC0 |> u"m"))
			@inbounds push!(dL, ustrip.(dL0 |> u"m"))
			@inbounds push!(dP, ustrip.(dP0 |> u"m"))
			@inbounds push!(dT, ustrip.(dT0 |> u"m"))
			@inbounds push!(dA, ustrip.(dA0 |> u"m"))
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

		# extrap_full = extrapolate(scale(interpolate(A, BSpline(Linear())), xs), Line())

		c2z = interpolate(dC[idxC], z[idxC], SteffenMonotonicInterpolation())
		l2z = interpolate(dL[idxL], z[idxL], SteffenMonotonicInterpolation())
		p2z = interpolate(dP[idxP], z[idxP], SteffenMonotonicInterpolation())
		t2z = interpolate(dT[idxT], z[idxT], SteffenMonotonicInterpolation())
		a2z = interpolate(dA[idxA], z[idxA], SteffenMonotonicInterpolation())

		return new{typeof(cosmo), typeof(h)}(h, ΩΛ, Ωm, Ωr, Ωk, wEOSΛ, cosmo, c2z, l2z, p2z, t2z, a2z)
	end
end

function CosmologicalModel(h::Real, Ωm::Real; z::Union{Nothing, Vector{Z}} = nothing) where {Z <: Real}
	cosmo = cosmology(h = h, OmegaM = Ωm)
	return CosmologicalModel{typeof(cosmo), typeof(h)}(cosmo; z = z)
end

function CosmologicalModel(h::Real, Ωm::Real, Ωk::Real; z::Union{Nothing, Vector{Z}} = nothing) where {Z <: Real}
	cosmo = cosmology(h = h, OmegaM = Ωm, OmegaK = Ωk)
	return CosmologicalModel{typeof(cosmo), typeof(h)}(cosmo; z = z)
end

function CosmologicalModel(h::Real, Ωm::Real, Ωk::Real, Ωr::Real; z::Union{Nothing, Vector{Z}} = nothing) where {Z <: Real}
	cosmo = cosmology(h = h, OmegaM = Ωm, OmegaK = Ωk, Omegar = Ωr)
	return CosmologicalModel{typeof(cosmo), typeof(h)}(cosmo; z = z)
end



# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Latest Planck's cosmology.
"""
const CosmologyPlanck(; z::Union{Nothing, Vector{Z}} = nothing) where {Z <: Real} = CosmologicalModel{typeof(cosmology()), Float64}(cosmology(); z = z)


# ----------------------------------------------------------------------------------------------- #
# 
# Overload Base functions. 
# The returned type is essentially the data type of the underlying cosmology from `Cosmology.jl` (Float64, by default).
Base.eltype(cosmo::CosmologicalModel) = typeof(cosmo._cosmology.h)



# ----------------------------------------------------------------------------------------------- #