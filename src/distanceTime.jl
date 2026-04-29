# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Abstract supertype for distance measurements.  
Sub-types include: 
- `DistanceComoving`
- `DistanceLightTravel`
- `DistanceAngularDiameter`
- `DistanceComovingTransverse`
- `DistanceLuminosity`
"""
abstract type AbstractDistanceMeasure end


# ----------------------------------------------------------------------------------------------- #
# 
@doc """
Abstract supertype for time measurements.
Sub-types include: 
- `TimeLookback`,
- `TimeConformal`.
"""
abstract type AbstractTimeMeasure end


# ----------------------------------------------------------------------------------------------- #
#
@doc """
	generateDistanceType(distanceTypeSym, info)

Generate a distance measure type named `Distance\$(distanceTypeSym)` together with all
standard constructors and conversions to/from `Redshift` and `ScaleFactor`.

# Input
- `distanceTypeSym` [`Symbol`]: suffix of the type name, e.g. `Comoving` → `DistanceComoving`
- `info` [`String`]: description text inserted into the type docstring
"""
macro generateDistanceType(distanceTypeSym, info)
	distanceType_str = string(distanceTypeSym)
	label = Symbol(lowercase(distanceType_str[1]) * distanceType_str[2:end])
	name = Symbol("Distance", distanceType_str)
	doc = """
		Convenient object to help with distance measures conversions.
		$(string(info))
		For more information see:
			"Distance measures in cosmology"
			D. Hogg
			arXiv:astro-ph/9905116

		# Input
		- `cosmology` [`CosmologicalModel`]: the cosmological model to be used as reference
		- `d` [`Length{T}`]: the distance
		"""
	return esc(quote
		mutable struct $(name){T <: Real} <: AbstractDistanceMeasure
			cosmology::CosmologicalModel
			value::Length{T}

			$(name){T}(cosmology::CosmologicalModel{C, T}, distance::Length) where {C, T} = begin
				val = T(ustrip(uconvert(u"Mpc", distance))) * u"Mpc"
				return new{T}(cosmology, val)
			end
		end
		@doc $(doc) $(name)

		$(name)(cosmology::CosmologicalModel{C, T}, distance::Length) where {C, T} = $(name){T}(cosmology, distance)
		$(name)(cosmology::CosmologicalModel{C, T}, distance::Real) where {C, T} = $(name){T}(cosmology, distance * u"Mpc")
		$(name)(cosmology::CosmologicalModel{C, T}, z::Redshift, z0::Redshift) where {C, T} = $(name)(cosmology, getfield(cosmology.fromRedshift, $(QuoteNode(label)))(z.value, z0.value))
		$(name)(cosmology::CosmologicalModel{C, T}, z::Redshift) where {C, T} = $(name)(cosmology, Redshift(T(0.)), z)
		$(name)(cosmology::CosmologicalModel{C, T}, a::ScaleFactor, a0::ScaleFactor) where {C, T} = $(name){T}(cosmology, Redshift(a), Redshift(a0))
		$(name)(cosmology::CosmologicalModel{C, T}, a::ScaleFactor) where {C, T} = $(name)(cosmology, ScaleFactor(a), ScaleFactor(T(1.)))
		$(name)(distance::Union{Length, Real}) = $(name)(defaultCosmology(), distance)
		$(name)(z::Redshift, z0::Redshift) = $(name)(defaultCosmology(), z, z0)
		$(name)(z::Redshift) = $(name)(defaultCosmology(), z)
		$(name)(a::ScaleFactor, a0::ScaleFactor) = $(name){T}(defaultCosmology(), a, a0)
		$(name)(a::ScaleFactor) = $(name)(defaultCosmology(), a)

		Redshift{Z}(distance::$(name)) where {Z} = Redshift{Z}(getfield(distance.cosmology.toRedshift, $(QuoteNode(label)))(distance.value))
		Redshift(distance::$(name)) = Redshift{eltype(distance)}(distance)

		ScaleFactor{A}(distance::$(name)) where {A} = Redshift{A}(distance) |> ScaleFactor
		ScaleFactor(distance::$(name)) = Redshift{eltype(distance)}(distance) |> ScaleFactor
	end)
end


@generateDistanceType(
	LightTravel,
	"""
	The light-travel distance corresponds to the time light from a given object would take to reach an observer.
	"""
)

@generateDistanceType(
	Comoving,
	"""
	The comoving distance is the distance between any two points in the reference frame of the Hubble flow.
	It is defined as:  
		d_c = R_H ∫ dz / E(z)
	between any two points z1 and z2.
	"""
)

@generateDistanceType(
	Luminosity,
	"""
	The luminosity distance is related to the flux (Φ) and the bolometric luminosity (L) as: d_L = sqrt(L / 4πΦ).
	"""
)


@generateDistanceType(
	AngularDiameter,
	"""
	The angular diameter distance is the ratio of an object's transverse length to its angular size.
	It relates to the transverse comoving distance in the following way:
	d_a = d_m / (1 + z)
	"""
)

@generateDistanceType(
	ComovingTransverse,
	"""
	For two objects at the same redshift separated by a given angle, the transverse comoving distance depends on the curvature:
	d_m = d_c	if  Ωk=0
	d_m = R_H sinh(sqrt(|Ωk|) d_c / R_H) / sqrt(|Ωk|) 	otherwise
	"""
)


# ----------------------------------------------------------------------------------------------- #
#
@doc """
	generateTimeType(timeTypeSym, info)

Generate a time measure type named `Time\$(timeTypeSym)` together with all standard
constructors and conversions to/from `Redshift` and `ScaleFactor`.

# Input
- `timeTypeSym` [`Symbol`]: suffix of the type name, e.g. `Lookback` → `TimeLookback`
- `info` [`String`]: description text inserted into the type docstring
"""
macro generateTimeType(timeTypeSym, info)
	timeType_str = string(timeTypeSym)
	label = Symbol(lowercase(timeType_str[1]) * timeType_str[2:end])
	name = Symbol("Time", timeType_str)
	doc = """
		Convenient object to help with time measures conversions.
		$(string(info))
		For more information see:
			"Distance measures in cosmology"
			D. Hogg
			arXiv:astro-ph/9905116

		# Input
		- `cosmology` [`CosmologicalModel`]: the cosmological model to be used as reference
		- `t` [`Time{T}`]: the time
		"""
	return esc(quote
		mutable struct $(name){T <: Real} <: AbstractTimeMeasure
			cosmology::CosmologicalModel
			value::Time{T}

			$(name){T}(cosmology::CosmologicalModel{C, T}, time::Time) where {C, T} = begin
				val = T(ustrip(uconvert(u"yr", time))) * u"yr"
				return new{T}(cosmology, val)
			end
		end
		@doc $(doc) $(name)

		$(name)(cosmology::CosmologicalModel{C, T}, time::Time) where {C, T} = $(name){T}(cosmology, time)
		$(name)(cosmology::CosmologicalModel{C, T}, time::Real) where {C, T} = $(name){T}(cosmology, time * u"yr")
		$(name)(cosmology::CosmologicalModel{C, T}, z::Redshift, z0::Redshift) where {C, T} = $(name)(cosmology, getfield(cosmology.fromRedshift, $(QuoteNode(label)))(z.value, z0.value))
		$(name)(cosmology::CosmologicalModel{C, T}, z::Redshift) where {C, T} = $(name)(cosmology, z, Redshift(T(0)))
		$(name)(cosmology::CosmologicalModel{C, T}, a::ScaleFactor, a0::ScaleFactor) where {C, T} = $(name)(cosmology, Redshift(a), Redshift(a0))
		$(name)(cosmology::CosmologicalModel{C, T}, a::ScaleFactor) where {C, T} = $(name)(cosmology, ScaleFactor(a), ScaleFactor(T(1.)))
		$(name)(time::Union{Time, Real}) = $(name)(defaultCosmology(), time)
		$(name)(z::Redshift, z0::Redshift) = $(name)(defaultCosmology(), z, z0)
		$(name)(z::Redshift) = $(name)(defaultCosmology(), z)
		$(name)(a::ScaleFactor, a0::ScaleFactor) = $(name){T}(defaultCosmology(), a, a0)
		$(name)(a::ScaleFactor) = $(name)(defaultCosmology(), a)

		Redshift{Z}(time::$(name)) where {Z} = Redshift{Z}(getfield(time.cosmology.toRedshift, $(QuoteNode(label)))(time.value))
		Redshift(time::$(name)) = Redshift{eltype(time)}(time)

		ScaleFactor{A}(time::$(name)) where {A} = Redshift{A}(time) |> ScaleFactor
		ScaleFactor(time::$(name)) = Redshift{eltype(time)}(time) |> ScaleFactor
	end)
end


@generateTimeType(
	Lookback,
	"""
	The lookback time corresponds to the time difference between the universe's age at two redshifts.
		t_l = t_H ∫ dz / E(z) / (1 + z)
	between any two points z1 and z2.
	"""
)


@generateTimeType(
	Conformal,
	"""
	The conformal time is the time in the frame of the Hubble flow.
	It is defined as:
		t_c = ∫ dt / a(t)
	between any two points z1 and z2.
	"""
)



# ----------------------------------------------------------------------------------------------- #
#                                  implementation of accessors                                    #
# ----------------------------------------------------------------------------------------------- #

for measure ∈ ("DistanceLightTravel", "DistanceComoving", "DistanceLuminosity", "DistanceAngularDiameter", "DistanceComovingTransverse", "TimeLookback", "TimeConformal")
	name = Symbol("$(measure)")
	@eval begin
		Base.eltype(::$(name){T}) where {T} = T
	end
end

# ----------------------------------------------------------------------------------------------- #
#                                    conversions and pipe                                         #
# ----------------------------------------------------------------------------------------------- #

for distanceType1 ∈ ("LightTravel", "Comoving", "Luminosity", "AngularDiameter", "ComovingTransverse")
	for distanceType2 ∈ ("LightTravel", "Comoving", "Luminosity", "AngularDiameter", "ComovingTransverse")
		d1 = Symbol("Distance$(distanceType1)")
		d2 = Symbol("Distance$(distanceType2)")

		if d1 ≠ d2
			@eval begin
				@doc """
				Constructor for `$($(d1)){D}` from `$($(d2))`.
				It essentially convert from `$($(d2))` to `$($(d1))`.
				It ultimately enables conversion implicit conversions and the usage of the operator `|>`.
				"""
				$(d1){D}(distance::($(d2))) where {D} = $(d1)(convert(D, distance.cosmology), Redshift{D}(distance))
				$(d1)(distance::($(d2))) = $(d1)(distance.cosmology, Redshift(distance))

				@doc """
					convert(DistanceType, distance)

				Type conversion from `$($(d2))` to `$($(d1))`.
				This function performs no type checks.
				It assumes that the underlying `CosmologicalModel` is the same for both.
				These conversions assume distances with respect to present time (z=0).
				"""
				Base.convert(::Type{$(d1){D1}}, d::$(d2)) where {D1} = $(d1)(convert(D1, d.cosmology), Redshift{D1}(d))
				Base.convert(::Type{$(d1)}, d::$(d2)) = $(d1)(d.cosmology, Redshift(d))
			end
		end	
	end
end

for timeType1 ∈ ("Lookback", "Conformal")
	for timeType2 ∈ ("Lookback", "Conformal")
		t1 = Symbol("Time$(timeType1)")
		t2 = Symbol("Time$(timeType2)")
		if t1 ≠ t2
			@eval begin

				@doc """
					convert(TimeType, time)

				Type conversion from `$($(t2))` to `$($(t1))`.
				It ultimately enables conversion implicit conversions and the usage of the operator `|>`.
				"""
				$(t1){T}(time::$(t2)) where {T} = $(t1)(convert(T, time.cosmology), Redshift{T}(time))
				$(t1)(time::($(t2))) = $(t1)(time.cosmology, Redshift(time))

				@doc """
					convert(TimeType, time)

				Type conversion from `$($(t2))` to `$($(t1))`.
				This function performs no type checks.
				It assumes that the underlying `CosmologicalModel` is the same for both.
				These conversions assume times with respect to present time (z=0).
				"""
				Base.convert(::Type{$(t1){T1}}, t::$(t2){T2}) where {T1, T2} = $(t1){promote_type(T1, T2)}(convert(T1, t.cosmology), Redshift(t))
				Base.convert(::Type{$(t1)}, t::$(t2){T2}) where {T2} = $(t2) |> $(t1)

			end
		end
	end
end



for measureType ∈ ("DistanceLightTravel", "DistanceComoving", "DistanceLuminosity", "DistanceAngularDiameter", "DistanceComovingTransverse", "TimeConformal", "TimeLookback")
	measure = Symbol("$(measureType)")
	@eval begin
		@doc """
			convert(Redshift, distance)

		Type conversion from `$($(measure))` to `Redshift`.
		These conversions assume distances/times with respect to present time (z=0).
		"""
		Base.convert(::Type{Redshift{Z}}, s::$(measure)) where {Z} = s |> Redshift{Z}
		Base.convert(::Type{Redshift}, s::$(measure)) = s |> Redshift
		

		@doc """
			convert(ScaleFactor, distance)

		Type conversion from `$($(measure))` to `ScaleFactor`.
		These conversions assume distances/times with respect to present time (z=0).
		"""
		Base.convert(::Type{ScaleFactor{A}}, s::$(measure)) where {A} = s |> ScaleFactor{A}
		Base.convert(::Type{ScaleFactor}, s::$(measure)) = s |> ScaleFactor


		@doc """
			convert(DistanceComoving, redshift, cosmology)
			convert(DistanceAngularDiameter, redshift, cosmology)
			convert(DistanceLightTravel, redshift, cosmology)
			convert(DistanceComovingTransverse, redshift, cosmology)
			convert(DistanceLuminosity, redshift, cosmology)

		Type conversion from `Redshift` to `$($(measure))`.
		These conversions assume distances/times with respect to present time (z=0).
		"""
		Base.convert(::Type{$(measure){S}}, z::Redshift, cosmology::CosmologicalModel) where {S} = $(measure)(convert(S, cosmology), convert(S, z))
		Base.convert(::Type{$(measure)}, z::Redshift, cosmology::CosmologicalModel) = convert($(measure){eltype(cosmology)}, z, cosmology)

		@doc """
			convert(DistanceComoving, scaleFactor, cosmology)
			convert(DistanceAngularDiameter, scaleFactor, cosmology)
			convert(DistanceLightTravel, scaleFactor, cosmology)
			convert(DistanceComovingTransverse, scaleFactor, cosmology)
			convert(DistanceLuminosity, scaleFactor, cosmology)

		Type conversion from `ScaleFactor` to `$($(measure))`.
		These conversions assume distances/times with respect to present time (z=0).
		"""
		Base.convert(::Type{$(measure){S}}, a::ScaleFactor, cosmology::CosmologicalModel) where {S} = $(measure)(convert(S, cosmology), a)
		Base.convert(::Type{$(measure)}, a::ScaleFactor, cosmology::CosmologicalModel) = convert($(measure){eltype(cosmology)}, a, cosmology)
	end
end

for measureType ∈ ("DistanceLightTravel", "DistanceComoving", "DistanceLuminosity", "DistanceAngularDiameter", "DistanceComovingTransverse", "TimeConformal", "TimeLookback")
	measure = Symbol("$(measureType)")
	@eval begin
		@doc """
			convert(DistanceType, distance)

		Type conversions for `$($(measure))`.
		This enables conversions between different distance measures (`DistanceLightTravel`, `DistanceComoving`, etc.) and the usage of the operator `|>`.
		"""
		Base.convert(::Type{$(measure){T}}, s::$(measure)) where {T} = $(measure){T}(convert(T, s.cosmology), T(s.value))
		Base.convert(::Type{T}, s::$(measure)) where {T <: Real} = $(measure){T}(convert(T, s.cosmology), T(s.value))

		Base.promote_rule(::Type{$(measure){T}}, ::Type{$(measure){U}}) where {T, U} = promote_type(T, U)
	end
end


# ----------------------------------------------------------------------------------------------- #
