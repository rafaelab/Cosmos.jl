using Cosmos, Unitful

cosmo = CosmologyPlanck()
time = TimeLookback(cosmo, Redshift(2.0))
conformal = TimeConformal(cosmo, time |> Redshift)

println("Time conversions")
println("  Lookback (yr): $(time.value)")
println("  Lookback (Gyr): $(uconvert(u"Gyr", time.value))")
println("  Redshift        : $(Redshift(time).value)")
println("  Conformal (yr)  : $(conformal.value)")
