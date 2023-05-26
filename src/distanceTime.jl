# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Abstract supertypes for distance and time measurements.
"""
abstract type AbstractDistanceMeasure end
abstract type AbstractTimeMeasure end


# ----------------------------------------------------------------------------------------------- #
# 
# Meta-generation of distance types.
for distanceType in ("LightTravel", "Comoving", "Luminosity", "AngularDiameter", "ComovingTransverse")
	label = Symbol(lowercase(distanceType[1]) * distanceType[2 : end])
	name = Symbol("Distance$(distanceType)")
	info = ""

	if distanceType == "LightTravel"
		info *= "The light-travel distance corresponds to the time light from a given object would take to reach an observer.\n"
	elseif distanceType == "Comoving"
		info *= "The comoving distance is the distance between any two points in the reference frame of the Hubble flow.\n"
		info *= "This is defined as: \n"
		info *= "  d_c = R_H ∫ dz / E(z) \n"
		info *= "between any two points z1 and z2.\n"
	elseif distanceType == "Luminosity"
		info *= "The luminosity distance is related to the flux (Φ) and the bolometric luminosity (L) as:\n"
		info *= "  d_L = sqrt(L / 4πΦ) .\n"
	elseif distanceType == "ComovingTransverse"
		info *= "For two objects at the same redshift separated by a given angle, the transverse comoving distance depends on the curvature:\n"
		info *= "  d_m = d_c	if  Ωk=0\n"
		info *= "  d_m = R_H sinh(sqrt(|Ωk|) d_c / R_H) / sqrt(|Ωk|) 	otherwise\n"
	elseif distanceType == "AngularDiameter"
		info *= "The angular diameter distance is the ratio of an object's transverse length to its angular size.\n"
		info *= "It relates to the transverse comoving distance in the following way:\n"
		info *= "  d_a = d_m / (1 + z)\n"
	end
	
	@eval begin
		@doc """
		Convenient object to help with distance measures conversions.
		$($(info))
		For more information see:
			"Distance measures in cosmology".
			D. Hogg
			arXiv:astro-ph/9905116

		# Input
		. `cosmology::CosmologicalModel`: the cosmological model to be used as reference \\
		. `d::Unitful.Length{T}`: the distance \\
		"""
		mutable struct $(name){T <: Real} <: AbstractDistanceMeasure
			cosmology::CosmologicalModel
			value::Unitful.Length{T}
			function $(name){T}(cosmology::CosmologicalModel{C, T}, distance::Unitful.Length) where {C, T}
				return new{T}(cosmology, distance |> u"Mpc")
			end
		end

		$(name)(cosmology::CosmologicalModel{C, T}, distance::Unitful.Length) where {C, T} = $(name){T}(cosmology, distance)

		$(name)(cosmology::CosmologicalModel{C, T}, distance::Real) where {C, T} = $(name){T}(cosmology, distance * u"Mpc")

		$(name)(cosmology::CosmologicalModel{C, T}, z::Redshift, z0::Redshift) where {C, T} = $(name)(cosmology, cosmology.fromRedshift[$(QuoteNode(label))](z0.value, z.value))

		$(name)(cosmology::CosmologicalModel{C, T}, z::Redshift) where {C, T} = $(name)(cosmology, z, Redshift(T(0.)))

		$(name)(cosmology::CosmologicalModel{C, T}, a::ScaleFactor, a0::ScaleFactor) where {C, T} = $(name){T}(cosmology, Redshift(a), Redshift(a0))

		$(name)(cosmology::CosmologicalModel{C, T}, a::ScaleFactor) where {C, T} = $(name)(cosmology, ScaleFactor(a), ScaleFactor(T(1.)))

		Redshift{Z}(distance::$(name)) where {Z} = Redshift{Z}(distance.cosmology.toRedshift[$(QuoteNode(label))](distance.value))

		Redshift(distance::$(name)) = Redshift{eltype(distance)}(distance)

		ScaleFactor{A}(distance::$(name)) where {A} = Redshift{A}(distance) |> ScaleFactor

		ScaleFactor(distance::$(name)) = Redshift{eltype(distance)}(distance) |> ScaleFactor
	end
end

# ----------------------------------------------------------------------------------------------- #
# 
# Meta-generation of time types.
for timeType in ("Lookback", "Conformal")
	label = Symbol(lowercase(timeType[1]) * timeType[2 : end])
	name = Symbol("Time$(timeType)")
	info = ""

	if timeType == "Lookback"
		info *= "The lookback time distance corresponds to the time difference between the universe's age at t2 and t1.\n"
		info *= "  t_l = t_H ∫ dz / E(z) / (1 + z) \n"
	elseif timeType == "Conformal"
		info *= "This is the time in the frame of the Hubble flow.\n"
		info *= "It is defined as: \n"
		info *= "  t_c = ∫ dt / a(t) \n"
		info *= "between any two points z1 and z2.\n"
	end
	
	@eval begin
		@doc """
		Convenient object to help with time measures conversions.
		$($(info))
		For more information see:
			"Distance measures in cosmology".
			D. Hogg
			arXiv:astro-ph/9905116

		# Input
		. `cosmology::CosmologicalModel`: the cosmological model to be used as reference \\
		. `t::Unitful.Time{T}`: the time \\
		"""
		mutable struct $(name){T <: Real} <: AbstractTimeMeasure
			cosmology::CosmologicalModel
			value::Unitful.Time{T}
			function $(name){T}(cosmology::CosmologicalModel{C, T}, time::Unitful.Time{U}) where {C, T, U}
				return new{T}(cosmology, time |> u"yr")
			end
		end

		$(name)(cosmology::CosmologicalModel{C, T}, time::Unitful.Length) where {C, T} = $(name){T}(cosmology, time)

		$(name)(cosmology::CosmologicalModel{C, T}, time::Real) where {C, T} = $(name){T}(cosmology, time * u"yr")

		$(name)(cosmology::CosmologicalModel{C, T}, z::Redshift, z0::Redshift) where {C, T} = $(name)(cosmology, cosmology.fromRedshift[$(QuoteNode(label))](z0.value, z.value))

		$(name)(cosmology::CosmologicalModel{C, T}, z::Redshift) where {C, T} = $(name)(cosmology, z, Redshift(T(0)))

		$(name)(cosmology::CosmologicalModel{C, T}, a::ScaleFactor, a0::ScaleFactor) where {C, T} = $(name)(cosmology, Redshift(a), Redshift(a0))

		$(name)(cosmology::CosmologicalModel{C, T}, a::ScaleFactor) where {C, T} = $(name)(cosmology, ScaleFactor(a), ScaleFactor(T(1.)))

		Redshift{Z}(time::$(name)) where {Z} = Redshift{Z}(time.cosmology.toRedshift[$(QuoteNode(label))](time.value))

		Redshift(time::$(name)) = Redshift{eltype(time)}(time)

		ScaleFactor{A}(time::$(name)) where {A} = Redshift{A}(time) |> ScaleFactor

		ScaleFactor(time::$(name)) = Redshift{eltype(time)}(time) |> ScaleFactor

	end
end

# ----------------------------------------------------------------------------------------------- #
#
# Implementation of `Base.eltype`.
for measure in ("DistanceLightTravel", "DistanceComoving", "DistanceLuminosity", "DistanceAngularDiameter", "DistanceComovingTransverse", "TimeLookback", "TimeConformal")
	name = Symbol("$(measure)")
	@eval begin
		Base.eltype(s::$(name){T}) where {T} = T
	end
end

# ----------------------------------------------------------------------------------------------- #
#
# Implement conversion functions between distance measures using `|>`.
for distanceType1 in ("LightTravel", "Comoving", "Luminosity", "AngularDiameter", "ComovingTransverse")
	for distanceType2 in ("LightTravel", "Comoving", "Luminosity", "AngularDiameter", "ComovingTransverse")
		d1 = Symbol("Distance$(distanceType1)")
		d2 = Symbol("Distance$(distanceType2)")
		if d1 ≠ d2
			@eval begin
				@doc """
				Type conversion from `$($(d2))` to `$($(d1))`.
				It ultimately enables conversion implicit conversions and the usage of the operator `|>`.
				"""
				$(d1)(distance::($(d2))) = $(d1)(distance.cosmology, Redshift(distance))
			end
		end
	end
end

# ----------------------------------------------------------------------------------------------- #
#
# Implement conversion functions between time measures using `|>`.
for timeType1 in ("Lookback", "Conformal")
	for timeType2 in ("Lookback", "Conformal")
		t1 = Symbol("Time$(timeType1)")
		t2 = Symbol("Time$(timeType2)")
		if t1 ≠ t2
			@eval begin
				@doc """
				Type conversion from `$($(t2))` to `$($(t1))`.
				It ultimately enables conversion implicit conversions and the usage of the operator `|>`.
				"""
				$(t1)(time::$(t2)) = $(t1)(time.cosmology, Redshift(time))
			end
		end
	end
end

# ----------------------------------------------------------------------------------------------- #
#
# Implement conversion functions between distance measures using `Base.convert`.
for distanceType1 in ("LightTravel", "Comoving", "Luminosity", "AngularDiameter", "ComovingTransverse")
	for distanceType2 in ("LightTravel", "Comoving", "Luminosity", "AngularDiameter", "ComovingTransverse")
		d1 = Symbol("Distance$(distanceType1)")
		d2 = Symbol("Distance$(distanceType2)")
		if d1 ≠ d2
			@eval begin
				@doc """
				Type conversion from `$($(d2))` to `$($(d1))`.
				This function performs no type checks.
				It assumes that the underlying `CosmologicalModel` is the same for both.
				These conversions assume distances with respect to present time (z=0).
				"""
				Base.convert(::Type{$(d1){D1}}, d::$(d2){D2}) where {D1, D2} = $(d1){promote_type(D1, D2)}(convert(D1, d.cosmology), Redshift(d))
				Base.convert(::Type{$(d1)}, d::$(d2){D2}) where {D2} = $(d2) |> $(d1)
			end
		end
	end
end

# ----------------------------------------------------------------------------------------------- #
#
# Implement conversion functions between time measures using `Base.convert`.
for timeType1 in ("Lookback", "Conformal")
	for timeType2 in ("Lookback", "Conformal")
		t1 = Symbol("Time$(timeType1)")
		t2 = Symbol("Time$(timeType2)")
		if t1 ≠ t2
			@eval begin
				@doc """
				Type conversion from `$($(t2))` to `$($(t1))`.
				This function performs no type checks.
				It assumes that the underlying `CosmologicalModel` is the same for both.
				These conversions assume times with respect to present time (z=0).
				"""
				Base.convert(::Type{$(t1){T1}}, t::$(t2){T2}) where {T1, T2} = $(t1){promote_type(T1, T2)}(convert(T1, d.cosmology), Redshift(t))
				Base.convert(::Type{$(t1)}, t::$(t2){T2}) where {T2} = $(t2) |> $(t1)
			end
		end
	end
end

# # ----------------------------------------------------------------------------------------------- #
# #
# @doc """
# Conversion from `TimeLookback` to `DistanceLightTravel` using Julia's built-in function `Base.convert`.
# """
# Base.convert(::Type{DistanceLightTravel}, time::TimeLookback) = (time.value * SpeedOfLightInVacuum) |> u"Mpc"

# # ----------------------------------------------------------------------------------------------- #
# #
# @doc """
# Conversion from `TimeConformal` to `DistanceComoving` using Julia's built-in function `Base.convert`.
# """
# Base.convert(::Type{DistanceComoving}, time::TimeConformal) = (time.value * SpeedOfLightInVacuum) |> u"Mpc"

# # ----------------------------------------------------------------------------------------------- #
# #
# @doc """
# Conversion from `DistanceLightTravel` to `TimeLookback` using Julia's built-in function `Base.convert`.
# """
# Base.convert(::Type{TimeLookback}, distance::DistanceLightTravel) = (distance.value / SpeedOfLightInVacuum) |> u"yr"

# # ----------------------------------------------------------------------------------------------- #
# #
# @doc """
# Conversion from `DistanceLightTravel` to `TimeLookback` using Julia's built-in function `Base.convert`.
# """
# Base.convert(::Type{TimeConformal}, distance::DistanceComoving) = (distance.value / SpeedOfLightInVacuum) |> u"yr"

# ----------------------------------------------------------------------------------------------- #

# DistanceComoving(cosmol::CosmologicalModel, z1::Redshift, z2::Redshift) = DistanceComoving{eltype(cosmol)}(cosmol.cosmology, comoving_radial_dist(cosmol.cosmology, z2.z, z1.z))
# DistanceComoving(cosmol::CosmologicalModel, z::Redshift) = DistanceComoving(cosmol, z, Redshift(0.))
	
# DistanceLuminosity(cosmol::CosmologicalModel, z::Redshift) = DistanceLuministy{eltype(cosmol)}(luminosity_dist(cosmol.cosmology, z.z))
# DistanceLuminosity(cosmol::CosmologicalModel, z1::Redshift, z2::Redshift) = DistanceLuminosity(cosmol, z1) - DistanceLuminosity(cosmol, z2)




# ----------------------------------------------------------------------------------------------- #
#

# Base.convert(::Type{Redshift{Z}}, d::DistanceLightTravel{D}) where {Z <: Real, D <: Real} = Redshift{Z}(d.cosmology._lightTravelDistance2Redshift(ustrip(d.d |> u"m")))
# Base.convert(::Type{Redshift{Z}}, d::DistanceLightTravel{D}) where {Z <: Real, D <: Real} = Redshift{Z}(d.cosmology._lightTravelDistance2Redshift(ustrip(d.d |> u"m")))

# redshiftToComovingDistance(cosmo::CosmologicalModel, z0::Real, z1::Real) = comoving_radial_dist(cosmo._cosmology, z0, z1) 
# redshiftToComovingDistance(cosmo::CosmologicalModel, z::Real) = comoving_radial_dist(cosmo._cosmology, zero(eltype(cosmo)), z)


# Redshift(d1::Light)



# ----------------------------------------------------------------------------------------------- #