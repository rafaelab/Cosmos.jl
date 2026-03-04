using Cosmos
using Unitful



cosmo = CosmologicalModel(0.68, 0.3, 0.01, 8e-5; Ωb = 0.048, Nν = 3.1)

println("Custom model summary")
println(cosmo)
println("Matter density (1/m^3):", computeMatterDensity(cosmo))
println("Baryon density (1/m^3):", computeBaryonDensity(cosmo))
