# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	@show(io, cosmology)
	print(io, cosmology)
	println(io, cosmology)

Display information of `CosmologicalModel` types.

# Input
. `io`: `IO`-type objects with standard output \\
. `cosmology`: `CosmologicalModel`-type object \\
"""
function Base.show(io::IO, cosmology::CosmologicalModel)
	printstyled(io, "$(typeof(cosmology)) = "; bold = true)
	if isflat(cosmology)
		print(io, "  . shape: flat\n")
	elseif isopen(cosmology)
		print(io, "  . shape: open\n")
		print(io, @sprintf("  . Ωk = %4.3e\n", cosmology.Ωk))
	elseif isclosed(cosmology)
		print(io, "  . shape: closed\n")
		print(io, @sprintf("  . Ωk = %4.3e\n", cosmology.Ωk))
	end
	print(io, @sprintf("  . ΩΛ = %4.3e\n", cosmology.ΩΛ))
	print(io, @sprintf("  . Ωm = %4.3e\n", cosmology.Ωm))
	print(io, @sprintf("  . Ωr = %4.3e\n", cosmology.Ωr))
	print(io, @sprintf("  . Ωb = %4.3e\n", cosmology.Ωb))
	print(io, @sprintf("  . h = %4.3f\n", cosmology.h))
	print(io, @sprintf("  . Tcmb = %4.3e\n", cosmology.Tcmb))
	print(io, @sprintf("  . Nν = %4.3f\n", cosmology.Nν))
	print(io, @sprintf("  . wEOSΛ(a) = %3.2f + %3.2f (1 - a)\n", cosmology.wEOSΛ...))
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	@show(io, z)
	print(io, z)
	println(io, z)

Display information of `Redshift` types.

# Input
. `io`: `IO`-type objects with standard output \\
. `z`: `Redshift`-type object \\
"""
function Base.show(io::IO, z::Redshift{Z}) where {Z}
	print(io, "z = $(z.value)")
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	@show(io, a)
	print(io, a)
	println(io, a)

Display information about `ScaleFactor` types.

# Input
. `io`: `IO`-type objects with standard output \\
. `a`: `ScaleFactor`-type object \\
"""
function Base.show(io::IO, a::ScaleFactor{A}) where {A}
	print(io, "a = $(a.value)")
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	@show(io, d)
	print(io, d)
	println(io, d)

Display information about `AbstractDistanceMeasure` types.

# Input
. `io`: `IO`-type objects with standard output \\
. `d`: `AbstractDistanceMeasure`-type object \\
"""
function Base.show(io::IO, d::AbstractDistanceMeasure)
	print(io, "d = $(d.value)")
end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	@show(io, t)
	print(io, t)
	println(io, t)

Display information about `AbstractTimeMeasure` types.

# Input
. `io`: `IO`-type objects with standard output \\
. `t`: `AbstractTimeMeasure`-type object \\
"""
function Base.show(io::IO, t::AbstractTimeMeasure)
	print(io, "t = $(t.value)")
end



# ----------------------------------------------------------------------------------------------- #
# 