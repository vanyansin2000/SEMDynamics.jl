@testset "halo periodic orbits" begin
    @test all(
        length(reference.seed) == 6
        for reference in values(SEMDynamics.PeriodicOrbits.ORBIT_REFERENCE_DATA)
    )
    references = Dict{Tuple{Symbol,Symbol},PeriodicOrbit}()

    for branch in (:northern, :southern), lp in (:L1, :L2)
        orbit = generate_halo(; branch, lp)
        references[(branch, lp)] = orbit

        @test orbit isa AbstractPeriodicOrbit
        @test length(orbit.x0) == 6
        @test orbit.x0[[2, 4, 6]] ≈ zeros(3) atol=1e-14
        @test !iszero(orbit.x0[3])
        @test isfinite(orbit.C)
        @test orbit.sol(orbit.P) ≈ orbit.x0 atol=1e-8
    end

    southern_nrho = generate_nrho_9_2(branch=:southern)
    northern_nrho = generate_nrho_9_2(branch=:northern)
    expected_nrho_period = 4pi / (9abs(Bcr4bp_Aux().EMRot.ws))
    direct_nrho = generate_halo(branch=:southern, lp=:L2, P=expected_nrho_period)
    short_period_halo = generate_halo(branch=:southern, lp=:L2, P=4pi / 9)
    @test southern_nrho.P ≈ northern_nrho.P ≈ expected_nrho_period
    @test southern_nrho.x0 ≈ direct_nrho.x0 atol=1e-9
    @test southern_nrho.x0[3] > 1e-3
    @test northern_nrho.x0[3] < -1e-3
    @test southern_nrho.sol(southern_nrho.P) ≈ southern_nrho.x0 atol=2e-7
    @test northern_nrho.sol(northern_nrho.P) ≈ northern_nrho.x0 atol=2e-7
    @test short_period_halo.P ≈ 4pi / 9
    @test short_period_halo.sol(short_period_halo.P) ≈ short_period_halo.x0 atol=2e-7

    corrected_dro = orbit_shooting(cr3bp_eqm!, generate_DRO().x0, pi)
    @test length(corrected_dro) == 6
    @test corrected_dro[3] == 0.0

    @test_throws ArgumentError generate_halo(branch=:eastern)
    @test_throws ArgumentError generate_halo(lp=:L3)
    @test_throws ArgumentError generate_halo(seed=[1.0, 0.1, 2.0])
    @test_throws ArgumentError orbit_shooting(cr3bp_eqm!, [1.0, 2.0], pi)
end
