# ----------------------------------------------------------------------------------------------- #
#
@doc """
	UnitInterpolation{I, U}

Callable struct pairing a monotonic interpolation object with a Unitful unit for input
normalisation.
Calling with a bare `Real` forwards the value directly to the interpolation (must already be
expressed in the correct SI unit).
Calling with a `Unitful.Quantity` first converts to `unit`, strips, then forwards.

# Members
- `interp` [`I`]: the underlying monotonic interpolation object
- `unit` [`U`]: Unitful unit used to normalise `Quantity` inputs (e.g. `u"m"`, `u"s"`)
"""
struct UnitInterpolation{I, U}
	interp::I
	unit::U
end

@inline (w::UnitInterpolation)(x::Real) = w.interp(x)
@inline (w::UnitInterpolation)(x::Unitful.Quantity) = w.interp(ustrip(uconvert(w.unit, x)))

# ----------------------------------------------------------------------------------------------- #
#
function sortedInterpolation(xs::Vector{T}, ys::Vector{T}) where {T <: Real}
	p = sortperm(xs)
	return interpolate(xs[p], ys[p], SteffenMonotonicInterpolation())
end

# ----------------------------------------------------------------------------------------------- #
#
@doc """
	RedshiftConversion{C, T, F}

Callable struct encoding one redshift-to-observable conversion.
Storing the specific Cosmology.jl function as the concrete type parameter `F` avoids
heap-allocated closures and enables full inlining.

# Members
- `cosmo` [`C`]: the underlying `AbstractCosmology` object
- `f` [`F`]: a top-level function `(cosmo, z2, z1) -> observable`
"""
struct RedshiftConversion{T <: Real, C <: AbstractCosmology, F}
	cosmo::C
	f::F
end

RedshiftConversion{T}(cosmo::C, f::F) where {T <: Real, C <: AbstractCosmology, F} = RedshiftConversion{T, C, F}(cosmo, f)


@inline (w::RedshiftConversion{T, C, F})(z2::Real, z1::Real) where {T, C, F} = T(w.f(w.cosmo, z2, z1))

# ----------------------------------------------------------------------------------------------- #
#
@inline _redshiftToComovingDistance(cosmo, z2::Real, z1::Real) = Cosmology.comoving_radial_dist(cosmo, z2, z1)
@inline _redshiftToLuminosityDistance(cosmo, z2::Real, z1::Real) = Cosmology.luminosity_dist(cosmo, z2) - Cosmology.luminosity_dist(cosmo, z1)
@inline _redshiftToTransverseComovingDistance(cosmo, z2::Real, z1::Real) = Cosmology.transverse_comoving_dist(cosmo, z2, z1)
@inline _redshiftToAngularDiameterDistance(cosmo, z2::Real, z1::Real) = Cosmology.angular_diameter_dist(cosmo, z1, z2)
@inline _redshiftToLookbackTime(cosmo, z2::Real, z1::Real) = Cosmology.lookback_time(cosmo, z2) - Cosmology.lookback_time(cosmo, z1)
@inline _redshiftToConformalTime(cosmo, z2::Real, z1::Real) = Cosmology.comoving_radial_dist(cosmo, z2, z1) / SpeedOfLightInVacuum |> u"yr"
@inline _redshiftToLightTravelDistance(cosmo, z2::Real, z1::Real) = (Cosmology.lookback_time(cosmo, z2) - Cosmology.lookback_time(cosmo, z1)) * SpeedOfLightInVacuum |> u"Mpc"

# ----------------------------------------------------------------------------------------------- #
# 
function prepareRedshiftSamples(T::Type{<: Real})
	z = T[]
	append!(z, -exp10.(range(-3., 0.; length = 91)))
	append!(z, collect(range(-0.001, 0.001; length = 31)))
	append!(z, collect(range(0.001, 10.; length = 81)))
	append!(z, collect(range(1., 2.; length = 21)))
	append!(z, exp10.(range(2., 4.; length = 41)))
	unique!(z)
	return z
end

# ----------------------------------------------------------------------------------------------- #
# 
function conversionsFromRedshift(cosmo::C, ::Type{T}) where {C <: AbstractCosmology, T <: Real}
	return (
		comoving = RedshiftConversion{T}(cosmo, _redshiftToComovingDistance),
		lightTravel = RedshiftConversion{T}(cosmo, _redshiftToLightTravelDistance),
		luminosity = RedshiftConversion{T}(cosmo, _redshiftToLuminosityDistance),
		transverseComoving = RedshiftConversion{T}(cosmo, _redshiftToTransverseComovingDistance),
		angularDiameter = RedshiftConversion{T}(cosmo, _redshiftToAngularDiameterDistance),
		lookback = RedshiftConversion{T}(cosmo, _redshiftToLookbackTime),
		conformal = RedshiftConversion{T}(cosmo, _redshiftToConformalTime),
	)
end

# ----------------------------------------------------------------------------------------------- #
# 
function conversionsToRedshift(cosmo::AbstractCosmology, z::AbstractVector, ::Type{T}) where {T <: Real}
	z = convert(Vector{T}, z)
	filter!(>(T(-1)), z)
	sort!(z)

	n = length(z)
	dC = Vector{T}(undef, n)
	dL = Vector{T}(undef, n)
	dP = Vector{T}(undef, n)
	dT = Vector{T}(undef, n)
	dA = Vector{T}(undef, n)
	tC = Vector{T}(undef, n)
	tL = Vector{T}(undef, n)

	Threads.@threads for i ∈ eachindex(z)
		dC0 = Cosmology.comoving_radial_dist(cosmo, z[i])
		dL0 = Cosmology.luminosity_dist(cosmo, z[i])
		dT0 = Cosmology.comoving_transverse_dist(cosmo, z[i])
		dA0 = Cosmology.angular_diameter_dist(cosmo, z[i])
		tL0 = Cosmology.lookback_time(cosmo, z[i])
		@inbounds dC[i] = ustrip(uconvert(u"m", dC0))
		@inbounds dL[i] = ustrip(uconvert(u"m", dL0))
		@inbounds dP[i] = ustrip(uconvert(u"m", tL0 * SpeedOfLightInVacuum))
		@inbounds dT[i] = ustrip(uconvert(u"m", dT0))
		@inbounds dA[i] = ustrip(uconvert(u"m", dA0))
		@inbounds tL[i] = ustrip(uconvert(u"s", tL0))
		@inbounds tC[i] = ustrip(uconvert(u"s", dC0 / SpeedOfLightInVacuum))
	end

	return (
		comoving = UnitInterpolation(sortedInterpolation(dC, z), u"m"),
		lightTravel = UnitInterpolation(sortedInterpolation(dP, z), u"m"),
		luminosity = UnitInterpolation(sortedInterpolation(dL, z), u"m"),
		transverseComoving = UnitInterpolation(sortedInterpolation(dT, z), u"m"),
		angularDiameter = UnitInterpolation(sortedInterpolation(dA, z), u"m"),
		lookback = UnitInterpolation(sortedInterpolation(tL, z), u"s"),
		conformal = UnitInterpolation(sortedInterpolation(tC, z), u"s"),
	)
end

# ----------------------------------------------------------------------------------------------- #
# 