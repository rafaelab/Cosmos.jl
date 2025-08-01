# Comparison with Ned Wright's CosmoCalc:
#   https://www.astro.ucla.edu/~wright/CosmoCalc.html
# Set h = 0.69, Omega_m = 0.29, Omega_lambda = 0.71
@testset "Comparison with CosmoCal - FlatLCDM" begin
	cosmo = CosmologyPlanck()

	# parameters for comparison
	z0 = Redshift(0.0)
	z1 = Redshift(0.1)
	z2 = Redshift(1.0)
	tolerance = 0.1 * u"Mpc"
	tolerance_t = 0.01 * u"Gyr"
	tolerance3 = 0.1 * u"Gpc ^ 3"
	dl1 = (1.321e9 * u"yr" |> u"s") * SpeedOfLightInVacuum |> u"Mpc"
	dl2 = (7.868e9 * u"yr" |> u"s") * SpeedOfLightInVacuum |> u"Mpc"
	dc1 = 424.8 * u"Mpc"
	dc2 = 3371.6 * u"Mpc"
	dl1 = 467.3 * u"Mpc"
	dl2 = 6743.1 * u"Mpc"
	da1 = 386.2 * u"Mpc"
	da2 = 1685.8 * u"Mpc"
	cv1 = 0.321 * u"Gpc ^ 3"
	cv2 = 160.543 * u"Gpc ^ 3"
	age0 = 13.787 * u"Gyr"
	age1 = 12.466 * u"Gyr"
	age2 = 5.918 * u"Gyr"

	@test DistanceLightTravel(cosmo, z1).value - dl1 < tolerance
	@test DistanceLightTravel(cosmo, z2).value - dl2 < tolerance
	@test DistanceComoving(cosmo, z1).value - dc1 < tolerance
	@test DistanceComoving(cosmo, z2).value - dc2 < tolerance
	@test DistanceLuminosity(cosmo, z1).value - dl1 < tolerance
	@test DistanceLuminosity(cosmo, z2).value - dl2 < tolerance
	@test DistanceAngularDiameter(cosmo, z1).value - dl1 < tolerance
	@test DistanceAngularDiameter(cosmo, z2).value - dl2 < tolerance
	# @test comovingVolume(cosmo, z1) - cv1 < tolerance3
	# @test comovingVolume(cosmo, z2) - cv2 < tolerance3
	# @test ageOfUniverse(cosmo, z0) - age0 < tolerance_t
	# @test ageOfUniverse(cosmo, z1) - age1 < tolerance_t
	# @test ageOfUniverse(cosmo, z2) - age2 < tolerance_t 
end

# Comparison with Ned Wright's CosmoCalc:
#   https://www.astro.ucla.edu/~wright/CosmoCalc.html
# Set h = 0.69, Omega_m = 0.29, Omega_lambda = 0.65, Omega_k = 0.06
@testset "Comparison with CosmoCal - OpenLCDM" begin
	cosmo = CosmologicalModel(0.69, 0.29, 0.06)

	# parameters for comparison
	z0 = 0.0
	z1 = 0.1
	z2 = 1.0
	tolerance = 0.1 * u"Mpc"
	tolerance_t = 0.01 * u"Gyr"
	tolerance3 = 0.1 * u"Gpc ^ 3"
	dc1 = 423.6 * u"Mpc"
	dc2 = 3311.9 * u"Mpc"
	dl1 = 466.9 * u"Mpc"
	dl2 = 7742. * u"Mpc"
	da1 = 385.1 * u"Mpc"
	da2 = 1665.6 * u"Mpc"
	cv1 = 0.318 * u"Gpc ^ 3"
	cv2 = 153.223 * u"Gpc ^ 3"
	age0 = 13.527 * u"Gyr"
	age1 = 12.209 * u"Gyr"
	age2 = 5.784 * u"Gyr"

	@test DistanceLightTravel(cosmo, z1).value - dl1 < tolerance
	@test DistanceLightTravel(cosmo, z2).value - dl2 < tolerance
	@test DistanceComoving(cosmo, z1).value - dc1 < tolerance
	@test DistanceComoving(cosmo, z2).value - dc2 < tolerance
	@test DistanceLuminosity(cosmo, z1).value - dl1 < tolerance
	@test DistanceLuminosity(cosmo, z2).value - dl2 < tolerance
	@test DistanceAngularDiameter(cosmo, z1).value - dl1 < tolerance
	@test DistanceAngularDiameter(cosmo, z2).value - dl2 < tolerance
	# @test comovingVolume(cosmo, z1) - cv1 < tolerance3
	# @test comovingVolume(cosmo, z2) - cv2 < tolerance3
	# @test ageOfUniverse(cosmo, z0) - age0 < tolerance_t
	# @test ageOfUniverse(cosmo, z1) - age1 < tolerance_t
	# @test ageOfUniverse(cosmo, z2) - age2 < tolerance_t 
end
