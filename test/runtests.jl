using Test
using SEMDynamics
using CairoMakie
using DifferentialEquations
using LinearAlgebra

include("visualization.jl")
include("periodic_orbits.jl")
include("dynamic_events.jl")

@testset "SEMDynamics.jl" begin
    @testset "parameters" begin
        aux = Bcr4bp_Aux()
        @test 0 < aux.EMRot.μ < 1
        @test aux.dim.EMRot_l > 0
        @test aux.EMRot.r_p1 > aux.EMRot.r_p2 > 0
    end

    @testset "CR3BP equations" begin
        aux = Bcr4bp_Aux()
        state2d = [0.8, 0.1, 0.02, -0.03]
        derivative2d = similar(state2d)
        cr3bp_eqm!(derivative2d, state2d, aux, 0.0)

        state3d = [state2d[1], state2d[2], 0.0, state2d[3], state2d[4], 0.0]
        derivative3d = similar(state3d)
        cr3bp_eqm!(derivative3d, state3d, aux, 0.0)

        @test all(isfinite, derivative2d)
        @test derivative3d ≈ [derivative2d[1], derivative2d[2], 0.0,
                              derivative2d[3], derivative2d[4], 0.0]

        state20 = [state2d; vec(Matrix{Float64}(I, 4, 4))]
        derivative20 = similar(state20)
        cr3bp_eqm!(derivative20, state20, aux, 0.0)
        @test derivative20[1:4] ≈ derivative2d
        @test all(isfinite, derivative20[5:end])

        state42 = [state3d; vec(Matrix{Float64}(I, 6, 6))]
        derivative42 = similar(state42)
        cr3bp_eqm!(derivative42, state42, aux, 0.0)
        @test derivative42[1:6] ≈ derivative3d
        @test all(isfinite, derivative42[7:end])
        @test_throws ErrorException cr3bp_eqm!(zeros(5), zeros(5), aux, 0.0)
    end

    @testset "BCR4BP equations" begin
        aux = Bcr4bp_Aux()
        state = [0.8, 0.1, 0.02, -0.03]
        derivative = similar(state)
        bcr4bp_eqm!(derivative, state, aux, 0.25)

        @test all(isfinite, derivative)
        @test derivative[1:2] == state[3:4]

        state20 = [state; vec(Matrix{Float64}(I, 4, 4))]
        derivative20 = similar(state20)
        bcr4bp_eqm!(derivative20, state20, aux, 0.25)
        @test derivative20[1:4] ≈ derivative
        @test all(isfinite, derivative20[5:end])

        state3d = [state[1], state[2], 0.0, state[3], state[4], 0.0]
        derivative3d = similar(state3d)
        bcr4bp_eqm!(derivative3d, state3d, aux, 0.25)
        state42 = [state3d; vec(Matrix{Float64}(I, 6, 6))]
        derivative42 = similar(state42)
        bcr4bp_eqm!(derivative42, state42, aux, 0.25)
        @test derivative42[1:6] ≈ derivative3d
        @test all(isfinite, derivative42[7:end])
        @test_throws ArgumentError bcr4bp_eqm!(zeros(5), zeros(5), aux, 0.0)
    end

    @testset "coordinate transformations" begin
        μ = Bcr4bp_Aux().EMRot.μ
        state = [0.2, -0.1, 0.03, 0.04]
        time = 0.37

        for center in (:p1, :p2)
            rotating = cr3bp_inertial_to_rotating(μ, time, state; center)
            recovered = cr3bp_rotating_to_inertial(μ, time, rotating; center)
            @test recovered ≈ state atol=1e-14
        end
        @test_throws ArgumentError cr3bp_inertial_to_rotating(μ, time, state; center=:sun)
    end

    @testset "energy and equilibrium points" begin
        aux = Bcr4bp_Aux()
        μ = aux.EMRot.μ
        state2d = [0.8, 0.1, 0.02, -0.03]
        state3d = [state2d[1], state2d[2], 0.0, state2d[3], state2d[4], 0.0]

        @test compute_jacobi(state2d, μ) ≈ compute_jacobi(state3d, μ)
        @test_throws ArgumentError compute_jacobi(zeros(5), μ)

        x_l1 = solve_L1_L2_x(μ; which=:L1)
        x_l2 = solve_L1_L2_x(μ; which=:L2)
        @test x_l1 < 1 - μ < x_l2

        for x in (x_l1, x_l2)
            derivative = zeros(4)
            cr3bp_eqm!(derivative, [x, 0.0, 0.0, 0.0], aux, 0.0)
            @test derivative ≈ zeros(4) atol=1e-10
        end
    end

    @testset "curve sampling" begin
        points = hcat(collect(0.0:4.0), zeros(5))
        @test uniform_sample_curve(points, 3; start_idx=1) == [1, 3, 5]
        @test_throws ArgumentError uniform_sample_curve(points, 0)
    end

    @testset "integration and events" begin
        parameters = ode_params(cr3bp_eqm!)
        final_state, times, states = integration(
            [0.8, 0.1, 0.02, -0.03],
            (0.0, 0.01),
            parameters;
            interp_num=5,
        )

        @test length(final_state) == 4
        @test length(times) == length(states) == 5
        @test all(isfinite, final_state)
        @test isempty(dynamic_events())
    end

    @testset "periodic orbit smoke test" begin
        orbit = generate_DRO(P=pi)
        @test orbit isa AbstractPeriodicOrbit
        @test orbit.P ≈ pi
        @test length(orbit.x0) == 6
        @test orbit.sol(orbit.P) ≈ orbit.x0 atol=1e-8
    end
end
