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
- `h` [`Real`]: dimensionless Hubble constant  
- `Ωm` [`Real`]: matter density  
- `Ωr` [`Real`]: radiation density   
- `Ωk` [`Real`]: curvature density  
- `ΩΛ` [`Real`]: dark energy density   
- `Ωb` [`Real`]: baryon density (set to -1 if unavailable)  
- `Nν` [`Real`]: number of effective neutrino species (defaults to 3)  
- `Tcmb` [`Real`]: CMB temperature at present time (defaults to 2.7255 K, following Planck)  
- `wEOSΛ` [`NTuple{2, Real}`]: tuple with parameters of the equation of state for dark energy: $w = w_0 + w_a (1 - a)$  
- `cosmology` [`AbstractCosmology``]: object from `Cosmology.jl`  
- `toRedshift` [`Dict{Symbol, Function}`]: functions to convert distance/time to redshift (`:comoving`, `:lightTravel`, `:angularDiameter`,`:transverseComoving`, `:luminosity`, `:lookback`, `:conformal`)  
- `fromRedshift::Dict{Symbol, Function}`: functions to convert distance/time from redshift (`:comoving`, `:lightTravel`, `:angularDiameter`, `:transverseComoving`, `:luminosity`, `:lookback`, `:conformal`)  
- `zArray` [`Vector{T}`]: array of values of redshift to build distance/time conversion functions; if nothing defaults to built-in values  


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
. Should this struct be immutable?
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
	toRedshift::ConversionTuple
	fromRedshift::ConversionTuple
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
function CosmologyPlanck(; T::Type = Float64)
	c0 = Cosmology.cosmology()
	cosmo = Cosmology.FlatLCDM{T}(T(c0.h), T(c0.Ω_Λ), T(c0.Ω_m), T(c0.Ω_r))
	return CosmologicalModel{typeof(cosmo), T}(cosmo; Ωb = T(0.023), Nν = T(3.04), Tcmb = T(2.7255))
end


# ----------------------------------------------------------------------------------------------- #
# 
setDefaultCosmology(CosmologyPlanck())


# ----------------------------------------------------------------------------------------------- #
# 
@doc """	
	isFlat(cosmology)

Determine whether a `CosmologicalModel` has a flat geometry.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isFlat(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractFlatCosmology


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isOpen(cosmology)	

Determine whether a `CosmologicalModel` has an open geometry.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isOpen(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractOpenCosmology


# ----------------------------------------------------------------------------------------------- #
# 
@doc """	
	isClosed(cosmology)

Determine whether a `CosmologicalModel` has a closed geometry.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isClosed(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractClosedCosmology


# ----------------------------------------------------------------------------------------------- #
# 
@doc """	
	isCold(cosmology)

Determine whether a `CosmologicalModel` is described by cold dark matter.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isCold(cosmology::CosmologicalModel) = begin 
	return 
		(cosmology.cosmology isa Cosmology.FlatLCDM) || 
		(cosmology.cosmology isa Cosmology.OpenLCDM) || 
		(cosmology.cosmology isa Cosmology.ClosedLCDM)
end


@doc """	
	isWarm(cosmology)

Determine whether a `CosmologicalModel` is described by warm dark matter.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isWarm(cosmology::CosmologicalModel) = begin 
	return 
		(cosmology.cosmology isa Cosmology.FlatWCDM) || 
		(cosmology.cosmology isa Cosmology.OpenWCDM) || 
		(cosmology.cosmology isa Cosmology.ClosedWCDM)
end


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
Base.:(==)(cosmol1::CosmologicalModel{C1, T1}, cosmol2::CosmologicalModel{C2, T2}) where {C1, C2, T1, T2} = begin
	(T1 == T2) && (C1 == C2) && (cosmol1.cosmology == cosmol2.cosmology) && (cosmol1.Ωb == cosmol2.Ωb) && (cosmol1.Nν == cosmol2.Nν) && (cosmol1.wEOSΛ == cosmol2.wEOSΛ)
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Object inequality comparison.
"""
Base.:(!=)(cosmol1::CosmologicalModel{C1, T1}, cosmol2::CosmologicalModel{C2, T2}) where {C1, C2, T1, T2} = ! (cosmol1 == cosmol2)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Type conversion.
"""
Base.convert(::Type{T}, cosmology::CosmologicalModel{C, U}) where {C, T <: Real, U} = begin 
	return CosmologicalModel{typeof(convert(T, cosmology.cosmology)), T}(convert(T, cosmology.cosmology); Ωb = cosmology.Ωb, Nν = cosmology.Nν, Tcmb = cosmology.Tcmb)
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Promotion rules.
"""
Base.promote_rule(::Type{CosmologicalModel{C, T1}}, ::Type{CosmologicalModel{C, T2}}) where {C, T1, T2} = CosmologicalModel{C, promote_type(T1, T2)}


# ----------------------------------------------------------------------------------------------- #
# 
