using Base.Threads
using Cosmology
import Cosmology: AbstractCosmology
using Cosmonstants: SpeedOfLightInVacuum
using Interpolations: SteffenMonotonicInterpolation, interpolate
using Unitful
using Unitful: u, ustrip, uconvert
using UnitfulAstro

function prepareRedshiftSamples(T::Type{<: Real})
	z = T[]
	append!(z, -10. .^ collect(range(-3., 0.; length = 91)))
	append!(z, collect(range(-0.001, 0.001; length = 31)))
	append!(z, collect(range(0.001, 10.; length = 81)))
	append!(z, collect(range(1., 2.; length = 21)))
	append!(z, 10 .^ collect(range(2., 4.; length = 41)))
	unique!(z)
	return z
end

function conversionsFromRedshift(cosmo::AbstractCosmology)
	T = eltype(cosmo)

	z2dl(z2::Real, z1::Real) = T(Cosmology.luminosity_dist(cosmo, z2) - Cosmology.luminosity_dist(cosmo, z1))
	z2dc(z2::Real, z1::Real) = T(Cosmology.comoving_radial_dist(cosmo, z2, z1))
	z2dt(z2::Real, z1::Real) = T(Cosmology.transverse_comoving_dist(cosmo, z2, z1))
	z2da(z2::Real, z1::Real) = T(Cosmology.angular_diameter_dist(cosmo, z1, z2))
	z2tl(z2::Real, z1::Real) = T(Cosmology.lookback_time(cosmo, z2) - Cosmology.lookback_time(cosmo, z1))
	z2tc(z2::Real, z1::Real) = T(z2dc(z2, z1) / SpeedOfLightInVacuum |> u"yr")
	z2dp(z2::Real, z1::Real) = T(z2tl(z2, z1) * SpeedOfLightInVacuum |> u"Mpc")

	return (
		comoving = z2dc,
		lightTravel = z2dp,
		luminosity = z2dl,
		transverseComoving = z2dt,
		angularDiameter = z2da,
		lookback = z2tl,
		conformal = z2tc
	)
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
		dC0 = Cosmology.comoving_radial_dist(cosmo, z[i])
		dL0 = Cosmology.luminosity_dist(cosmo, z[i])
		dT0 = Cosmology.comoving_transverse_dist(cosmo, z[i])
		dA0 = Cosmology.angular_diameter_dist(cosmo, z[i])
		tL0 = Cosmology.lookback_time(cosmo, z[i])
		dP0 = tL0 * SpeedOfLightInVacuum
		tC0 = dC0 / SpeedOfLightInVacuum
		@inbounds dC[i] = ustrip.(uconvert(u"m", dC0))
		@inbounds dL[i] = ustrip.(uconvert(u"m", dL0))
		@inbounds dP[i] = ustrip.(uconvert(u"m", dP0))
		@inbounds dT[i] = ustrip.(uconvert(u"m", dT0))
		@inbounds dA[i] = ustrip.(uconvert(u"m", dA0))
		@inbounds tL[i] = ustrip.(uconvert(u"s", tL0))
		@inbounds tC[i] = ustrip.(uconvert(u"s", tC0))
	end

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
	dc2z(x::Unitful.Length) = dc2z_(ustrip(uconvert(u"m", x)))
	dl2z(x::Unitful.Length) = dl2z_(ustrip(uconvert(u"m", x)))
	dp2z(x::Unitful.Length) = dp2z_(ustrip(uconvert(u"m", x)))
	dt2z(x::Unitful.Length) = dt2z_(ustrip(uconvert(u"m", x)))
	da2z(x::Unitful.Length) = da2z_(ustrip(uconvert(u"m", x)))
	tl2z(x::Unitful.Length) = tl2z_(ustrip(uconvert(u"m", x)))
	tc2z(x::Unitful.Length) = tc2z_(ustrip(uconvert(u"m", x)))

	return (
		comoving = dc2z,
		lightTravel = dp2z,
		luminosity = dl2z,
		transverseComoving = dt2z,
		angularDiameter = da2z,
		lookback = tl2z,
		conformal = tc2z
	)
end
