@doc """	isFlat(cosmology)

Determine whether a `CosmologicalModel` has a flat geometry.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isFlat(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractFlatCosmology

@doc """	isOpen(cosmology)

Determine whether a `CosmologicalModel` has an open geometry.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isOpen(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractOpenCosmology

@doc """	isClosed(cosmology)

Determine whether a `CosmologicalModel` has a closed geometry.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isClosed(cosmology::CosmologicalModel) = cosmology.cosmology isa Cosmology.AbstractClosedCosmology

@doc """	isCold(cosmology)

Determine whether a `CosmologicalModel` is described by cold dark matter.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isCold(cosmology::CosmologicalModel) = (cosmology.cosmology isa Cosmology.FlatLCDM) || (cosmology.cosmology isa Cosmology.OpenLCDM) || (cosmology.cosmology isa Cosmology.ClosedLCDM)

@doc """	isWarm(cosmology)

Determine whether a `CosmologicalModel` is described by warm dark matter.

# Input
- `cosmol::CosmologicalModel`: the cosmological model
"""
isWarm(cosmology::CosmologicalModel) = (cosmology.cosmology isa Cosmology.FlatWCDM) || (cosmology.cosmology isa Cosmology.OpenWCDM) || (cosmology.cosmology isa Cosmology.ClosedWCDM)

@doc """
Get type of values contained in `CosmologicalModel` object.
"""
Base.eltype(cosmology::CosmologicalModel) = typeof(cosmology.h)

@doc """
Object equality comparison.
"""
Base.:(==)(cosmol1::CosmologicalModel{C1, T1}, cosmol2::CosmologicalModel{C2, T2}) where {C1, C2, T1, T2} = begin
	(T1 == T2) && (C1 == C2) && (cosmol1.cosmology == cosmol2.cosmology) && (cosmol1.Ωb == cosmol2.Ωb) && (cosmol1.Nν == cosmol2.Nν) && (cosmol1.wEOSΛ == cosmol2.wEOSΛ)
end

Base.:(!=)(cosmol1::CosmologicalModel{C1, T1}, cosmol2::CosmologicalModel{C2, T2}) where {C1, C2, T1, T2} = ! (cosmol1 == cosmol2)

@doc """
Type conversion.
"""
Base.convert(::Type{T}, cosmology::CosmologicalModel{C, U}) where {C, T <: Real, U} = begin 
	return CosmologicalModel{typeof(convert(T, cosmology.cosmology)), T}(convert(T, cosmology.cosmology); Ωb = cosmology.Ωb, Nν = cosmology.Nν, Tcmb = cosmology.Tcmb)
end

@doc """
Promotion rules.
"""
Base.promote_rule(::Type{CosmologicalModel{C, T1}}, ::Type{CosmologicalModel{C, T2}}) where {C, T1, T2} = CosmologicalModel{C, promote_type(T1, T2)}
