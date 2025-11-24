# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	@show(io, cosmology)
	print(io, cosmology)
	println(io, cosmology)

Display information of `CosmologicalModel` types.


# Input
- `io` [`IO`]: `IO`-type objects with standard output 
- `cosmology` [`CosmologicalModel`]: the cosmological model object 
"""
function Base.show(io::IO, cosmology::CosmologicalModel)
	printstyled(io, "$(typeof(cosmology)) = \n"; bold = true)
	if isFlat(cosmology)
		print(io, "  . shape: flat\n")
	elseif isOpen(cosmology)
		print(io, "  . shape: open\n")
	elseif isClosed(cosmology)
		print(io, "  . shape: closed\n")
	end
	print(io, @sprintf("  . ΩΛ = %4.3e\n", cosmology.ΩΛ))
	print(io, @sprintf("  . Ωk = %4.3e\n", cosmology.Ωk))
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
- `io` [`IO`]: `IO`-type objects with standard output 
- `z` [`Redshift`]: the redshift 
"""
Base.show(io::IO, z::Redshift) = print(io, "z = $(z.value)")


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	@show(io, a)
	print(io, a)
	println(io, a)

Display information about `ScaleFactor` types.


# Input
- `io` [`IO`]: `IO`-type objects with standard output 
- `a` [`ScaleFactor`]: the scale factor 
"""
Base.show(io::IO, a::ScaleFactor) = print(io, "a = $(a.value)")


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	@show(io, d)
	print(io, d)
	println(io, d)

Display information about `AbstractDistanceMeasure` types.
Note that distances are cosmology-dependent, but this information is not displayed.


# Input
- `io` [`IO`]: `IO`-type objects with standard output 
- `d` [`AbstractDistanceMeasure`]: a distance object 
"""
Base.show(io::IO, d::AbstractDistanceMeasure) = print(io, "d = $(d.value)")


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	@show(io, t)
	print(io, t)
	println(io, t)

Display information about `AbstractTimeMeasure` types.
Note that times are cosmology-dependent, but this information is not displayed.


# Input
- `io` [`IO`]: `IO`-type objects with standard output 
- `t` [`AbstractTimeMeasure`]: the time measure 
"""
Base.show(io::IO, t::AbstractTimeMeasure) = print(io, "t = $(t.value)")



# ----------------------------------------------------------------------------------------------- #
# 