# ----------------------------------------------------------------------------------------------- #
# 
@doc """
	isLengthDimension(d)

Check if dimension provided is correct.
"""
@inline function isLengthDimension(d::Unitful.AbstractQuantity) 
	dimension(d) == Unitful.𝐋 || throw(DimensionMismatch("Dimension of provided quantity is not distance."))
	return true
end

# ----------------------------------------------------------------------------------------------- #
# 