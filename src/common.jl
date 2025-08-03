# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Useful alias for type definition.
To prevent conflicts, this is not exported.
"""
const Maybe{T} = Union{Nothing, T}


# ----------------------------------------------------------------------------------------------- #
# 
for cosmologyType ∈ ("FlatLCDM", "OpenLCDM", "ClosedLCDM", "FlatWCDM", "OpenWCDM", "ClosedWCDM")
	@eval begin
		@doc """
		Useful extensions of `Cosmology.jl`` to enable type conversion related to `$($(Symbol("$(cosmologyType)")))`.
		"""
		Base.convert(::Type{T}, cosmo::$(Symbol("$(cosmologyType)")){U}) where {T <: Real, U} = begin
			fields = fieldnames(typeof(cosmo))
			values = Tuple([convert(T, getfield(cosmo, field)) for field in fields])
			return $(Symbol("$(cosmologyType)")){T}(values...)
		end

		Base.convert(::Type{$(Symbol("$(cosmologyType)")){T}}, cosmo::$(Symbol("$(cosmologyType)")){U}) where {T, U} = convert(T, cosmo)

		@doc """
		Useful extensions of Cosmology.jl implementing promotion rules.
		"""
		Base.promote_rule(::Type{$(Symbol("$(cosmologyType)")){T}}, ::Type{$(Symbol("$(cosmologyType)")){U}}) where {T, U} = begin
			return $(Symbol("$(cosmologyType)")){promote_type(T, U)}
		end
	end
end

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Get the underlying data type of an `AbstractCosmology` object (from `Cosmology.jl`).
"""
Base.eltype(cosmo::AbstractCosmology) = typeof(cosmo.h)


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Define global variable to hold information about the default cosmological model.
"""
const defaultCosmologyRef = Ref{Any}(nothing)

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Set the value of the global default cosmological model.
This will enable faster function calls. 
For instances, instead of `d = DistanceComovingTransverse(1., cosmology)`, the second argument will become the default value.
"""
setDefaultCosmology(cosmology) = (defaultCosmologyRef[] = cosmology)
getDefaultCosmology() = defaultCosmologyRef[]


# ----------------------------------------------------------------------------------------------- #