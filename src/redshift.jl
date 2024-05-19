# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Convenient object to handle redshifts.
"""
struct Redshift{T <: Real}
	value::T
	function Redshift{T}(z::Real) where {T}
		z > -1. || throw(DomainError("Redshift cannot be less than or equal to -1."))
		return new{T}(z)
	end
end

Redshift(z::Real) = Redshift{typeof(z)}(z)

# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Convenient object to handle scale factors.
"""
struct ScaleFactor{T <: Real}
	value::T
	function ScaleFactor{T}(a::Real) where {T}
		a ≥ 0. || throw(DomainError("Scale factor cannot be negative or zero."))
		return new{T}(a)
	end
end

ScaleFactor(a::Real) = ScaleFactor{typeof(a)}(a)


# ----------------------------------------------------------------------------------------------- #
#
@doc """
	eltype(redshift)
	eltype(scaleFactor)

Get data type.
"""
Base.eltype(::Redshift{Z}) where {Z} = Z
Base.eltype(::ScaleFactor{A}) where {A} = A

# ----------------------------------------------------------------------------------------------- #
#
@doc """
Conversions for redshift and scale factors.
"""
Base.convert(::Type{Redshift{Z1}}, z::Redshift{Z2}) where {Z1, Z2} = Redshift{Z1}(Z1(z.value))
Base.convert(::Type{Redshift{Z}}, scaleFactor::ScaleFactor{A}) where {Z, A} = Redshift{Z}(1. / scaleFactor.value - 1.)
Base.convert(::Type{Redshift}, scaleFactor::ScaleFactor{A}) where {A} = convert(Redshift{A}, scaleFactor)
Base.convert(::Type{Z}, z::Redshift) where {Z <: Real} = Redshift{Z}(z.value)

Base.convert(::Type{ScaleFactor{A1}}, a::ScaleFactor{A2}) where {A1, A2} = ScaleFactor{A1}(A1(a.value))
Base.convert(::Type{ScaleFactor{A}}, redshift::Redshift{Z}) where {Z, A} = ScaleFactor{A}(1. / (1. + redshift.value))
Base.convert(::Type{ScaleFactor}, redshift::Redshift{Z}) where {Z} = convert(ScaleFactor{Z}, redshift)
Base.convert(::Type{A}, a::ScaleFactor) where {A <: Real} = ScaleFactor{A}(a.value)

Base.promote_rule(::Type{Redshift{Z1}}, ::Type{Redshift{Z2}}) where {Z1, Z2} = Redshift{promote_type(Z1, Z2)}
Base.promote_rule(::Type{ScaleFactor{A1}}, ::Type{ScaleFactor{A2}}) where {A1, A2} = ScaleFactor{promote_type(A1, A2)}

ScaleFactor(redshift::Redshift) = convert(ScaleFactor{eltype(redshift)}, redshift)
Redshift(scaleFactor::ScaleFactor) = convert(Redshift{eltype(scaleFactor)}, scaleFactor)


# ----------------------------------------------------------------------------------------------- #
#
@doc """
Operations with `Redshift`.
"""
Base.:(==)(z1::Redshift, z2::Redshift) = z1.value == z2.value
Base.:(!=)(z1::Redshift, z2::Redshift) = z1.value ≠ z2.value
Base.isequal(z1::Redshift, z2::Redshift) = z1.value === z2.value


# ----------------------------------------------------------------------------------------------- #
#
@doc """
Operations with `ScaleFactor`.
"""
Base.:(==)(a1::ScaleFactor, a2::ScaleFactor) = a1.value == a2.value
Base.:(!=)(a1::ScaleFactor, a2::ScaleFactor) = a1.value ≠ a2.value
Base.isequal(a1::ScaleFactor, a2::ScaleFactor) = a1.value === a2.value


# ----------------------------------------------------------------------------------------------- #
#