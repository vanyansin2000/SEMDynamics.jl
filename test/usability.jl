@testset "usability and failure contracts" begin
    aux = Bcr4bp_Aux()
    μ = aux.EMRot.μ
    planar = [0.8, 0.1, 0.02, -0.03]
    spatial = [0.8, 0.1, 0.07, 0.02, -0.03, 0.04]
    t = 0.37

    for center in (:p1, :p2)
        inertial = cr3bp_rotating_to_inertial(μ, t, spatial; center)
        @test cr3bp_inertial_to_rotating(μ, t, inertial; center) ≈ spatial atol=1e-14
        # At t=0 the inertial velocity includes ω × r about the selected primary.
        offset = center === :p1 ? μ : μ - 1
        expected = [spatial[1] + offset, spatial[2], spatial[3],
                    spatial[4] - spatial[2], spatial[5] + spatial[1] + offset, spatial[6]]
        @test cr3bp_rotating_to_inertial(μ, 0.0, spatial; center) ≈ expected
    end

    for state in (planar, spatial)
        n = length(state)
        extended = [state; vec(Matrix{Float64}(I, n, n))]
        inertial = cr3bp_rotating_to_inertial(μ, t, state; center=:p2)
        dimension = n ÷ 2
        expected_energy = sum(abs2, inertial[dimension+1:end]) / 2 - μ / norm(inertial[1:dimension])
        @test compute_ε2(state, t, μ) ≈ expected_energy
        @test compute_ε2(extended, t, μ) ≈ expected_energy
        @test compute_ε2(state, t + 0.5, μ) ≈ expected_energy
        for dynamics in (cr3bp_eqm!, bcr4bp_eqm!)
            du = similar(state)
            dynamics(du, state, aux, t)
            h = 1e-6
            finite_difference = (compute_ε2(state + h * du, t + h, μ) -
                                 compute_ε2(state - h * du, t - h, μ)) / (2h)
            rate = compute_ε2_dot(state, t, μ, aux; dynamics)
            @test rate ≈ finite_difference atol=1e-8 rtol=1e-6
            @test compute_ε2_dot(extended, t, μ, aux; dynamics) ≈ rate
        end
        θ = aux.EMRot.ws * t
        H, Hdot, Γ = SEMDynamics.Dynamics.Hamiltonian(state, θ, aux)
        @test all(isfinite, (H, Hdot, Γ))
        du = similar(state)
        bcr4bp_eqm!(du, state, aux, t)
        h = 1e-6
        finite_difference = (
            first(SEMDynamics.Dynamics.Hamiltonian(state + h * du, θ + h * aux.EMRot.ws, aux)) -
            first(SEMDynamics.Dynamics.Hamiltonian(state - h * du, θ - h * aux.EMRot.ws, aux))
        ) / (2h)
        @test Hdot ≈ finite_difference atol=1e-5 rtol=1e-5
        @test SEMDynamics.Dynamics.energy_zero_condition(extended, t, (; p=aux)) ≈ expected_energy
    end

    for callback in (cb_enter, cb_escape, cb_apse_p1, cb_apse_p2, cb_perilune, cb_apolune)
        @test callback(dynamic_events(); root_find=false).rootfind ==
              callback(dynamic_events(); rootfind=false).rootfind
    end

    function constant_velocity!(du, u, p, t)
        du[1] = u[3]
        du[2] = u[4]
        du[3] = du[4] = 0.0
        return nothing
    end
    for (primary_x, radius, callback) in ((-μ, aux.EMRot.r_p1, cb_p1collision),
                                        (1-μ, aux.EMRot.r_p2, cb_p2collision))
        events = dynamic_events()
        state = [primary_x + 1.55radius, 0.0, -radius, 0.0]
        parameters = ode_params(constant_velocity!, (; adaptive=false, dt=0.1), aux)
        @test integration(state, (0.0, 1.0), parameters; cb=callback(events)) ==
              (nothing, nothing, nothing)
        empty!(events)
        sol = integration(state, (0.0, 1.0), parameters; cb=callback(events), return_solution=true)
        @test sol.retcode == ReturnCode.Terminated
        @test sol.t[end] ≈ 0.6 atol=1e-12
        @test length(events) == 1
        @test events[1].state == sol.u[end]
    end

    # Energy plots use the ODE's own model and parameters, including the solar force.
    bcr_solution = solve(ODEProblem(bcr4bp_eqm!, spatial, (t, t + 0.01), aux), Vern7();
                         abstol=1e-12, reltol=1e-12)
    figure = Figure()
    @test length(plot_traj_ε2!(Axis(figure[1, 1]), bcr_solution,
                              [Event(:origin, t, copy(spatial))])) == 2

    parameters = ode_params(cr3bp_eqm!)
    sol = integration(planar, (0.0, 0.01), parameters; return_solution=true)
    @test sol.retcode == ReturnCode.Success
    # A solver failure remains inspectable in the opt-in native result.
    failed = integration(planar, (0.0, 10.0), parameters;
                         odeargs=(; maxiters=1), return_solution=true)
    @test !SEMDynamics.Dynamics.SciMLBase.successful_retcode(failed)

    orbit = generate_DRO()
    @test first(orbit.sol.t) == -orbit.P
    @test last(orbit.sol.t) == 2orbit.P
    @test orbit.sol(-orbit.P) ≈ orbit.x0 atol=1e-8
    @test orbit.sol(0.0) ≈ orbit.x0 atol=1e-8
    @test_throws ODESolveFailedError SEMDynamics.PeriodicOrbits.integrate_orbit(
        cr3bp_eqm!, orbit.x0, orbit.P; maxiters=1)
end
