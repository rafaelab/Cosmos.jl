@testset "Default cosmology and conversions" begin
	Cosmos.setDefaultCosmology(CosmologyPlanck())
	z = Redshift(0.5)
	dc = DistanceComoving(z)

	@test Cosmos.defaultCosmology() isa CosmologicalModel
	@test dc.cosmology === Cosmos.defaultCosmology()
	@test Redshift(dc).value ≈ z.value atol = 1e-5
end

@testset "Baryon density guard" begin
	cosmo = CosmologicalModel(0.69, 0.29)
	@test_throws ArgumentError computeBaryonDensity(cosmo)
end

@testset "Time conversion robustness" begin
	cosmo = CosmologyPlanck()
	tl = TimeLookback(cosmo, 1.0u"Gyr")
	@test uconvert(u"Gyr", tl.value) ≈ 1.0 * u"Gyr"

	tc = convert(TimeConformal{Float32}, TimeConformal(cosmo, Redshift(0.5)))
	@test eltype(tc) == Float32
end
