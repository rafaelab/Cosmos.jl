# ----------------------------------------------------------------------------------------------- #
#
@doc """
General struct to hold a `Cosmology.jl` object of type `AbstractCosmology`.
It exposes the relevant cosmological parameters and adds new ones.
The equation of state follows the Chevallier-Polarski-Linder parametrisation: \\
	"Accelerating Universes with Scaling Dark Matter" \\
	M. Chevallier and D. Polarski \\
	International Journal of Modern Physics D 10 (2001) 213. \\
	https://arxiv.org/abs/gr-qc/0009008 \\
	https://doi.org/10.1142/S0218271801000822 \\
	\\
	"Exploring the Expansion History of the Universe" \\
	E. Linder \\
	Physical Review Letters 90 (2003) 091301. \\
	https://doi.org/10.1103/PhysRevLett.90.091301 \\
	\\
The default constructors can be built using only the first 3 or 4 parameters.

# Members
. `h` [`Real`]: dimensionless Hubble constant \\
. `Ωm` [`Real`]: matter density \\
. `Ωr` [`Real`]: radiation density \\
. `Ωk` [`Real`]: curvature density \\
. `ΩΛ` [`Real`]: dark energy density \\
. `Ωb` [`Real`]: baryon density (set to -1 if unavailable) \\
. `Nν` [`Real`]: number of effective neutrino species (defaults to 3) \\
. `Tcmb` [`Real`]: CMB temperature at present time (defaults to 2.7255 K, following Planck) \\
. `wEOSΛ` [`NTuple{2, Real}`]: tuple with parameters of the equation of state for dark energy: `w = w_0 + w_a (1 - a)` \\
. `cosmology` [`AbstractCosmology``]: object from `Cosmology.jl` \\
. `toRedshift` [`Dict{Symbol, Function}`]: functions to convert distance/time to redshift (`:comoving`, `:lightTravel`, `:angularDiameter`,`:transverseComoving`, `:luminosity`, `:lookback`, `:conformal`) \\
. `fromRedshift::Dict{Symbol, Function}`: functions to convert distance/time from redshift (`:comoving`, `:lightTravel`, `:angularDiameter`, `:transverseComoving`, `:luminosity`, `:lookback`, `:conformal`) \\
. `zArray` [`Vector{T}`]: array of values of redshift to build distance/time conversion functions; if nothing defaults to built-in values \\

# Examples
```
	# define parameters
	Tcmb = 2.7255
	h = 0.69
	ΩΛ = 0.7099
	Ωk = 0.
	Ωm = 0.29
	Ωr = 1. - ΩΛ - Ωk - Ωm
	Nν = 3.04
	
	# some constructors
	cosmo1 = CosmologicalModel(Cosmology.FlatLCDM{Float64}(h, ΩΛ, Ωm, Ωr); Nν = Nν, Tcmb = Tcmb)
	cosmo2 = CosmologicalModel(h, Ωm; Tcmb = Tcmb,  Nν = Nν) # assumes Ωr = 0
	cosmo3 = CosmologicalModel(h, Ωm, Ωk; Tcmb = Tcmb,  Nν = Nν) # if geometry is not flat and Ωr = 0
	cosmo4 = CosmologicalModel(h, Ωm, Ωk, Ωr; Tcmb = Tcmb,  Nν = Nν) # includes radiation and non-flat geometry 
```

# To do
. Consider taking `Unitful` quantities. \\
. Should this struct be immutable?\\
"""
mutable struct CosmologicalModel{C <: AbstractCosmology, T <: Real}
	h::T
	ΩΛ::T
	Ωm::T
	Ωr::T
	Ωk::T
	Ωb::T
	Tcmb::T
	Nν::T
	wEOSΛ::SVector{2, T}
	cosmology::C
	toRedshift::Dict{Symbol, Function}
	fromRedshift::Dict{Symbol, Function}
	zArray::Vector{T}

	function CosmologicalModel{C, T}(cosmo::C; Ωb::Real = -1., Tcmb::Real = 2.7255, Nν::Real = 3., z::Maybe{AbstractVector} = nothing) where {C <: AbstractCosmology, T <: Real}
		h = T(cosmo.h)
		Ωr = T(cosmo.Ω_r)
		Ωm = T(cosmo.Ω_m)
		ΩΛ = T(cosmo.Ω_Λ)
		Ωb = T(Ωb)
		Ωk = (C isa Union{Cosmology.AbstractClosedCosmology, Cosmology.AbstractOpenCosmology}) ? T(cosmo.Ω_k) : zero(T)
		wEOSΛ = (C isa Union{Cosmology.FlatWCDM, Cosmology.OpenWCDM, Cosmology.ClosedWCDM}) ? SVector{2, T}(T(cosmo.w0), T(cosmo.wa)) : SVector{2, T}(-one(T), zero(T))
		Tcmb = T(Tcmb)
		Nν = T(Nν)

		if isnothing(z)
			z = prepareRedshiftSamples(T)
		end

		funcTo = conversionsToRedshift(cosmo, z)
		funcFrom = conversionsFromRedshift(cosmo)

		return new{C, T}(h, ΩΛ, Ωm, Ωr, Ωk, Ωb, Tcmb, Nν, wEOSΛ, convert(T, cosmo), funcTo, funcFrom, z)
	end
end

CosmologicalModel{T}(cosmo::AbstractCosmology; args...) where {T <: Real} = begin
	c = convert(T, cosmo)
	return CosmologicalModel{typeof(c), T}(c; args...)
end

CosmologicalModel(cosmo::AbstractCosmology; args...) = begin
	return CosmologicalModel{typeof(c), eltype(cosmo)}(cosmo; args...)
end

CosmologicalModel(h::Real, Ωm::Real; args...) = begin
	cosmo = Cosmology.cosmology(h = h, OmegaM = Ωm)
	return CosmologicalModel{typeof(cosmo), typeof(h)}(cosmo; args...)
end

CosmologicalModel(h::Real, Ωm::Real, Ωk::Real; args...) = begin
	cosmo = Cosmology.cosmology(h = h, OmegaM = Ωm, OmegaK = Ωk)
	return CosmologicalModel{typeof(cosmo), typeof(h)}(cosmo; args...)
end

CosmologicalModel(h::Real, Ωm::Real, Ωk::Real, Ωr::Real; args...) = begin
	cosmo = Cosmology.cosmology(h = h, OmegaM = Ωm, OmegaK = Ωk, OmegaR = Ωr)
	return CosmologicalModel{typeof(cosmo), typeof(h)}(cosmo; args...)
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Helper function (unexported) to iniatilise conversions of distance or time to/from redshift.

# To do
. Consider speeding up building this function using `@threads`
"""
function conversionsFromRedshift(cosmo::AbstractCosmology)
	T = eltype(cosmo)

	z2dl(z2::Real, z1::Real) = T(luminosity_dist(cosmo, z2) - luminosity_dist(cosmo, z1))
	z2dc(z2::Real, z1::Real) = T(comoving_radial_dist(cosmo, z2, z1))
	z2dt(z2::Real, z1::Real) = T(transverse_comoving_dist(cosmo, z2, z1))
	z2da(z2::Real, z1::Real) = T(angular_diameter_dist(cosmo, z1, z2))
	z2tl(z2::Real, z1::Real) = T(lookback_time(cosmo, z2) - lookback_time(cosmo, z1) )
	z2tc(z2::Real, z1::Real) = T(z2dc(z2, z1) / SpeedOfLightInVacuum |> u"yr")
	z2dp(z2::Real, z1::Real) = T(z2tl(z2, z1) * SpeedOfLightInVacuum |> u"Mpc")

	funcDict = Dict()
	funcDict[:comoving] = z2dc
	funcDict[:lightTravel] = z2dp
	funcDict[:luminosity] = z2dl
	funcDict[:transverseComoving] = z2dt
	funcDict[:angularDiameter] = z2da
	funcDict[:lookback] = z2tl
	funcDict[:conformal] = z2tc

	return funcDict
end

function conversionsToRedshift(cosmo::AbstractCosmology, z::AbstractVector)
	T = eltype(cosmo)
	z = convert(Vector{T}, z)
	z = z[z .> -1.0]
	sort!(z)

	# pre-allocate arrays
	dC = Vector{T}(undef, length(z))
	dL = Vector{T}(undef, length(z))
	dP = Vector{T}(undef, length(z))
	dT = Vector{T}(undef, length(z))
	dA = Vector{T}(undef, length(z))
	tC = Vector{T}(undef, length(z))
	tL = Vector{T}(undef, length(z))


	Threads.@threads for i ∈ eachindex(z)
		dC0 = comoving_radial_dist(cosmo, z[i])
		dL0 = luminosity_dist(cosmo, z[i]) 
		dT0 = comoving_transverse_dist(cosmo, z[i])
		dA0 = angular_diameter_dist(cosmo, z[i])
		tL0 = lookback_time(cosmo, z[i])
		dP0 = tL0 * SpeedOfLightInVacuum
		tC0 = dC0 / SpeedOfLightInVacuum
		@inbounds dC[i] = ustrip.(dC0 |> u"m")
		@inbounds dL[i] = ustrip.(dL0 |> u"m")
		@inbounds dP[i] = ustrip.(dP0 |> u"m")
		@inbounds dT[i] = ustrip.(dT0 |> u"m")
		@inbounds dA[i] = ustrip.(dA0 |> u"m")
		@inbounds tL[i] = ustrip.(tL0 |> u"s")
		@inbounds tC[i] = ustrip.(tC0 |> u"s")
	end

	# get indices that sort arrays, as required by Interpolations.jl
	idxDC = sortperm(dC)
	idxDL = sortperm(dL)
	idxDP = sortperm(dP)
	idxDT = sortperm(dT)
	idxDA = sortperm(dA)
	idxTC = sortperm(tC)
	idxTL = sortperm(tL)

	dc2z_ = interpolate(dC[idxDC], z[idxDC], SteffenMonotonicInterpolation())
	dl2z_ = interpolate(dL[idxDL], z[idxDL], SteffenMonotonicInterpolation())
	dp2z_ = interpolate(dP[idxDP], z[idxDP], SteffenMonotonicInterpolation())
	dt2z_ = interpolate(dT[idxDT], z[idxDT], SteffenMonotonicInterpolation())
	da2z_ = interpolate(dA[idxDA], z[idxDA], SteffenMonotonicInterpolation())
	tl2z_ = interpolate(tL[idxTL], z[idxTL], SteffenMonotonicInterpolation())
	tc2z_ = interpolate(tC[idxTC], z[idxTC], SteffenMonotonicInterpolation())

	dc2z(x::Real) = dc2z_(x)
	dl2z(x::Real) = dl2z_(x)
	dp2z(x::Real) = dp2z_(x)
	dt2z(x::Real) = dt2z_(x)
	da2z(x::Real) = da2z_(x)
	tl2z(x::Real) = tl2z_(x)
	tc2z(x::Real) = tc2z_(x)
	dc2z(x::Unitful.Length) = dc2z_(ustrip(x |> u"m"))
	dl2z(x::Unitful.Length) = dl2z_(ustrip(x |> u"m"))
	dp2z(x::Unitful.Length) = dp2z_(ustrip(x |> u"m"))
	dt2z(x::Unitful.Length) = dt2z_(ustrip(x |> u"m"))
	da2z(x::Unitful.Length) = da2z_(ustrip(x |> u"m"))
	tl2z(x::Unitful.Length) = tl2z_(ustrip(x |> u"m"))
	tc2z(x::Unitful.Length) = tc2z_(ustrip(x |> u"m"))

	funcDict = Dict()
	funcDict[:comoving] = dc2z
	funcDict[:lightTravel] = dp2z
	funcDict[:luminosity] = dl2z
	funcDict[:transverseComoving] = dt2z
	funcDict[:angularDiameter] = da2z
	funcDict[:lookback] = tl2z
	funcDict[:conformal] = tc2z

	return funcDict
end

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Latest Planck's cosmology.
This is based on:
	"Planck 2018 results. VI. Cosmological parameters"
	Planck Collaboration
	Astronomy and Astrophysics 641 (2020) A6.
	https://arxiv.org/abs/1807.06209
	https://doi.org/10.1051/0004-6361/201833910

# Input
. `z::Maybe{AbstractVector}`: vectors at which the redshifts will be sampled (for interpolating distance measures) \\
. `T::Type`: type of the data (defaults to `Float64`) \\
"""
function CosmologyPlanck(; T::Type = Float64)
	cosmo0 = Cosmology.cosmology()
	cosmo = Cosmology.FlatLCDM{T}(T(cosmo0.h), T(cosmo0.Ω_Λ), T(cosmo0.Ω_m), T(cosmo0.Ω_r))
	return CosmologicalModel{typeof(cosmo), T}(cosmo; Ωb = T(0.023), Nν = T(3.04), Tcmb = T(2.7255))
end

setDefaultCosmology(CosmologyPlanck())

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isFlat(cosmology)

Determine whether a `CosmologicalModel` has a flat geometry.

# Input
. `cosmol::CosmologicalModel`: the cosmological model
"""
isFlat(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractFlatCosmology


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isOpen(cosmology)

Determine whether a `CosmologicalModel` has an open geometry.

# Input
. `cosmol::CosmologicalModel`: the cosmological model
"""
isOpen(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractOpenCosmology


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isClosed(cosmology)

Determine whether a `CosmologicalModel` has a closed geometry.

# Input
. `cosmol::CosmologicalModel`: the cosmological model
"""
isClosed(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractClosedCosmology


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isCold(cosmology)

Determine whether a `CosmologicalModel` is described by cold dark matter.

# Input
. `cosmol::CosmologicalModel`: the cosmological model
"""
isCold(cosmology::CosmologicalModel) = (cosmology.cosmology isa Cosmology.FlatLCDM) || (cosmology.cosmology isa Cosmology.OpenLCDM) || (cosmology.cosmology isa Cosmology.ClosedLCDM)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isWarm(cosmology)

Determine whether a `CosmologicalModel` is described by warm dark matter.

# Input
. `cosmol::CosmologicalModel`: the cosmological model
"""
isWarm(cosmology::CosmologicalModel) = (cosmology.cosmology isa Cosmology.FlatWCDM) || (cosmology.cosmology isa Cosmology.OpenWCDM) || (cosmology.cosmology isa Cosmology.ClosedWCDM)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Get type of values contained in `CosmologicalModel` object.
"""
Base.eltype(cosmology::CosmologicalModel) = typeof(cosmology.h)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Object equality comparison.
"""
Base.:(==)(cosmol1::CosmologicalModel{C1, T1}, cosmol2::CosmologicalModel{C2, T2}) where {C1, C2, T1, T2} = (T1 == T2) && (C1 == C2) && (cosmol1.cosmology == cosmol2.cosmology) && (cosmol1.Ωb == cosmol2.Ωb) && (cosmol1.Nν == cosmol2.Nν) && (cosmol1.wEOSΛ .== cosmol1.wEOSΛ)

Base.:(!=)(cosmol1::CosmologicalModel{C1, T1}, cosmol2::CosmologicalModel{C2, T2}) where {C1, C2, T1, T2} = !(cosmol1 == cosmol2)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Type conversion.
"""
Base.convert(::Type{T}, cosmology::CosmologicalModel{C, U}) where {C, T <: Real, U} = CosmologicalModel{typeof(convert(T, cosmology.cosmology)), T}(convert(T, cosmology.cosmology); Ωb = cosmology.Ωb, Nν = cosmology.Nν, Tcmb = cosmology.Tcmb)

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Promotion rules.
"""
Base.promote_rule(::Type{CosmologicalModel{C, T1}}, ::Type{CosmologicalModel{C, T2}}) where {C, T1, T2} = CosmologicalModel{C, promote_type(T1, T2)}


# ----------------------------------------------------------------------------------------------- #