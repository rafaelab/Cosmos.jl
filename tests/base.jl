@testset "Check initialisation and property storage" begin
	cosmo1 = CosmologyPlanck()
	cosmo2 = CosmologicalModel(0.69, 0.29, 0., 0.06)

	tol1 = 1e-3
	tol2 = 1e-3
	
	@test cosmo1.h == 0.69
	@test cosmo1.Ωm == 0.29
	@test cosmo1.Ωk == 0.00
	@test cosmo1.Ωr - 8.778e-5 < tol1
	@test cosmo1.ΩΛ - (1. - cosmo1.Ωm - cosmo1.Ωr) < tol1
	@test cosmo2.h == 0.69
	@test cosmo2.Ωm == 0.29
	# @test cosmo2.Ωk == 0.06
	# @test cosmo2.Ωr == 0.00
	@test cosmo2.ΩΛ - (1. - cosmo2.Ωm - cosmo2.Ωr - cosmo2.Ωk) < tol2
end


@testset "Check geometry" begin
	planck = CosmologyPlanck()
	@test isFlat(planck)
	@test ! isOpen(planck)
	@test ! isClosed(planck)
end

@testset "Type conversion and promotion" begin
	planckF64 = CosmologyPlanck()
	planckF32 = convert(Float32, planckF64)
	@test eltype(planckF32.h) == Float32
end