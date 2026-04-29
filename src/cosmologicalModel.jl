# ----------------------------------------------------------------------------------------------- #
# 
@doc raw"""
General struct to hold a `Cosmology.jl` object of type `AbstractCosmology`.
It exposes the relevant cosmological parameters and adds new ones.
The equation of state follows the Chevallier-Polarski-Linder parametrisation:  
	"Accelerating Universes with Scaling Dark Matter"
	M. Chevallier and D. Polarski
	International Journal of Modern Physics D 10 (2001) 213.
	https://arxiv.org/abs/gr-qc/0009008
	https://doi.org/10.1142/S0218271801000822
	"Exploring the Expansion History of the Universe"
	E. Linder
	Physical Review Letters 90 (2003) 091301.
	https://doi.org/10.1103/PhysRevLett.90.091301  

The default constructors can be built using only the first 3 or 4 parameters.  


# Members
- `Ωb` [`T`]: baryon density parameter (set to -1 if unavailable)
- `Tcmb` [`T`]: CMB temperature at present time in Kelvin (defaults to 2.7255 K, following Planck)
- `Nν` [`T`]: number of effective neutrino species (defaults to 3)
- `wEOSΛ` [`SVector{2, T}`]: equation-of-state parameters for dark energy: $w(a) = w_0 + w_a (1 - a)$; set to $(-1, 0)$ for ΛCDM
- `cosmology` [`C`]: underlying `Cosmology.jl` object (concrete type, e.g. `FlatLCDM{T}`)
- `fromRedshift` [`FF`]: named tuple of concrete closures `(z2, z1) -> distance/time` (keys: `:comoving`, `:lightTravel`, `:luminosity`, `:transverseComoving`, `:angularDiameter`, `:lookback`, `:conformal`)
- `toRedshift` [`TF`]: named tuple of concrete closures `distance/time -> z` (same keys)

# Indirect members
- `h`: the dimensionless Hubble parameter, delegated to the underlying `cosmology` field
- `ΩΛ`: the dark energy density parameter, delegated to the underlying `cosmology` field
- `Ωm`: the matter density parameter, delegated to the underlying `cosmology` field
- `Ωr`: the radiation density parameter, delegated to the underlying `cosmology` field
- `Ωk`: the curvature density parameter, delegated to the underlying `cosmology` field (returns 0 for flat cosmologies)


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
. Consider taking `Unitful` quantities.
"""
struct CosmologicalModel{C <: AbstractCosmology, T <: Real, FF, TF}
	Ωb::T
	Tcmb::T
	Nν::T
	wEOSΛ::SVector{2, T}
	cosmology::C
	fromRedshift::FF
	toRedshift::TF

	CosmologicalModel{C, T}(cosmo::C; Ωb::Real = -1., Tcmb::Real = 2.7255, Nν::Real = 3., z::Maybe{AbstractVector} = nothing) where {C <: AbstractCosmology, T <: Real} = begin
		Ωb_ = T(Ωb)
		Tcmb_ = T(Tcmb)
		Nν_ = T(Nν)
		wEOSΛ = (C <: Union{Cosmology.FlatWCDM, Cosmology.OpenWCDM, Cosmology.ClosedWCDM}) ? SVector{2, T}(T(cosmo.w0), T(cosmo.wa)) : SVector{2, T}(-one(T), zero(T))
		zSamples = isnothing(z) ? prepareRedshiftSamples(T) : convert(Vector{T}, z)
		funcFrom = conversionsFromRedshift(cosmo, T)
		funcTo = conversionsToRedshift(cosmo, zSamples, T)
		return new{C, T, typeof(funcFrom), typeof(funcTo)}(Ωb_, Tcmb_, Nν_, wEOSΛ, convert(T, cosmo), funcFrom, funcTo)
	end
end

CosmologicalModel{T}(cosmo::AbstractCosmology; args...) where {T <: Real} = begin
	c = convert(T, cosmo)
	return CosmologicalModel{typeof(c), T}(c; args...)
end

CosmologicalModel(cosmo::AbstractCosmology; args...) = begin
	return CosmologicalModel{typeof(cosmo), eltype(cosmo)}(cosmo; args...)
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
	getproperty(cosmology, name)

Property access for `CosmologicalModel`.
The fields `h`, `ΩΛ`, `Ωm`, `Ωr`, `Ωk` are not stored directly; they are delegated to the underlying `cosmology` field to avoid data duplication.
Accessing them in a tight loop is fine; for absolute maximum performance use `cosmology.cosmology.h` etc. to bypass dispatch.
"""
function Base.getproperty(cosmology::CosmologicalModel{C, T}, name::Symbol) where {C, T}
	c = getfield(cosmology, :cosmology)
	if name === :h
		return c.h
	elseif name === :ΩΛ
		return c.Ω_Λ
	elseif name === :Ωm
		return c.Ω_m
	elseif name === :Ωr
		return c.Ω_r
	elseif name === :Ωk
		return C <: Union{Cosmology.AbstractClosedCosmology, Cosmology.AbstractOpenCosmology} ? c.Ω_k : zero(T)
	else
		return getfield(cosmology, name)
	end
end

@inline Base.propertynames(::CosmologicalModel) = (:Ωb, :Tcmb, :Nν, :wEOSΛ, :cosmology, :fromRedshift, :toRedshift, :h, :ΩΛ, :Ωm, :Ωr, :Ωk)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	CosmologyPlanck()

Convenience constructor for the Planck 2018 cosmological parameters (TT,TE,EE+lowE+lensing).

# Reference
"Planck 2018 results. VI. Cosmological parameters"
Planck Collaboration
Astronomy and Astrophysics 641 (2020) A6
arXiv:1807.06209
doi:10.1051/0004-6361/201833910
bibkey: planck2020a
"""
function CosmologyPlanck(; T::Type = Float64)
	c0 = Cosmology.cosmology()
	cosmo = Cosmology.FlatLCDM{T}(T(c0.h), T(c0.Ω_Λ), T(c0.Ω_m), T(c0.Ω_r))
	return CosmologicalModel{typeof(cosmo), T}(cosmo; Ωb = T(0.023), Nν = T(3.04), Tcmb = T(2.7255))
end



# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isFlat(cosmology)

Determine whether a `CosmologicalModel` has a flat geometry.

# Input
- `cosmology` [`CosmologicalModel`]: the cosmological model
"""
@inline isFlat(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractFlatCosmology


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isOpen(cosmology)	

Determine whether a `CosmologicalModel` has an open geometry.

# Input
- `cosmology` [`CosmologicalModel`]: the cosmological model
"""
@inline isOpen(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractOpenCosmology



# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isClosed(cosmology)

Determine whether a `CosmologicalModel` has a closed geometry.

# Input
- `cosmology` [`CosmologicalModel`]: the cosmological model
"""
@inline isClosed(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractClosedCosmology


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isCold(cosmology)

Determine whether a `CosmologicalModel` is described by cold dark matter.

# Input
- `cosmology` [`CosmologicalModel`]: the cosmological model
"""
@inline isCold(::CosmologicalModel{C}) where {C} = C <: Union{Cosmology.FlatLCDM, Cosmology.OpenLCDM, Cosmology.ClosedLCDM}


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isWarm(cosmology)

Determine whether a `CosmologicalModel` is described by warm dark matter.

# Input
- `cosmology` [`CosmologicalModel`]: the cosmological model
"""
@inline isWarm(cosmology::CosmologicalModel{C}) where {C} = C <: Union{Cosmology.FlatWCDM, Cosmology.OpenWCDM, Cosmology.ClosedWCDM}


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Get type of values contained in `CosmologicalModel` object.
"""
@inline Base.eltype(::CosmologicalModel{C, T}) where {C, T} = T


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Object equality comparison.
"""
function Base.:(==)(cosmol1::CosmologicalModel, cosmol2::CosmologicalModel)
	return getfield(cosmol1, :cosmology) == getfield(cosmol2, :cosmology) && cosmol1.Ωb == cosmol2.Ωb && cosmol1.Nν == cosmol2.Nν && cosmol1.wEOSΛ == cosmol2.wEOSΛ
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Type conversion.
"""
function Base.convert(::Type{T}, cosmology::CosmologicalModel{C, U}) where {C, T <: Real, U}
	return CosmologicalModel{typeof(convert(T, cosmology.cosmology)), T}(convert(T, cosmology.cosmology); Ωb = cosmology.Ωb, Nν = cosmology.Nν, Tcmb = cosmology.Tcmb)
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Promotion rules.
"""
function Base.promote_rule(::Type{CosmologicalModel{C, T1, FF1, TF1}}, ::Type{CosmologicalModel{C, T2, FF2, TF2}}) where {C, T1, T2, FF1, TF1, FF2, TF2}
	return CosmologicalModel{C, promote_type(T1, T2)}
end


# ----------------------------------------------------------------------------------------------- #
# 
