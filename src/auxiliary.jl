# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	prepareRedshiftSamples(T::Type{<: Real})

Prepare a set of redshift samples for cosmological calculations.
This function generates a range of redshift values, including negative, small, and large redshifts,
to cover a wide range of cosmological scenarios.
The redshift values are returned as a unique array of type `T`.
"""
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

# ----------------------------------------------------------------------------------------------- #
# 