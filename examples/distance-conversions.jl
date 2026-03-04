using Cosmos, Unitful

cosmo = CosmologyPlanck()
distance = DistanceComoving(cosmo, 600u"Mpc")
redshift = distance |> Redshift
scale = redshift |> ScaleFactor
luminosity = convert(DistanceLuminosity, distance)

println("Planck distance conversions for z=$(redshift.value)")
println("  Comoving   : $(distance.value)")
println("  Scale factor: $(scale.value)")
println("  Luminosity  : $(luminosity.value)")
