@testset "Cosmology functions" begin
    # Define a sample cosmological model
    Tcmb = 2.7255
    h = 0.69
    ΩΛ = 0.7099
    Ωk = 0.0
    Ωm = 0.29
    Ωr = 1.0 - ΩΛ - Ωk - Ωm
    Nν = 3.04
    cosmo = Cosmology.cosmology(; h = h, OmegaM = Ωm, OmegaK = Ωk, OmegaR = Ωr)
    model = CosmologicalModel(cosmo; Nν = Nν, Tcmb = Tcmb)

#     @testset "Hubble Parameter Functions" begin
#         @test dimensionlessHubbleParameter(model, 0.0) ≈ 1.0
#         @test hubbleParameter(model, 0.0) ≈ hubbleConstant(model)
#         @test hubbleDistance(model) > 0.0
#         @test hubbleTime(model) > 0.0
#     end

#     @testset "Scale Factor" begin
#         @test scaleFactor(model, 0.0) ≈ 1.0
#         @test scaleFactor(1.0) ≈ 0.5
#     end

#     @testset "Age of Universe" begin
#         @test ageOfUniverse(model) > 0.0
#         @test ageOfUniverse(model, 0.0) ≈ ageOfUniverse(model)
#     end

#     @testset "Comoving Volume and Element" begin
#         @test comovingVolume(model, 1.0) > 0.0
#         @test comovingVolumeElement(model, 1.0) > 0.0
#         @test comovingElement(model, 1.0) > 0.0
#     end

#     @testset "Density Calculations" begin
#         @test calculateDensityPhotons(model, Tcmb) > 0.0
#         @test calculateDensityNeutrinos(model, Nν, Tcmb) > 0.0

#     end
end