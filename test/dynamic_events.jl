@testset "dynamic event callbacks" begin
    aux = Bcr4bp_Aux()
    μ = aux.EMRot.μ
    sphere_radius = aux.EMRot.d2lim
    p2_x = 1 - μ

    function constant_velocity!(du, u, p, t)
        du[1] = u[4]
        du[2] = u[5]
        du[3] = u[6]
        du[4] = 0
        du[5] = 0
        du[6] = 0
        return nothing
    end

    enter_events = dynamic_events()
    enter_state = [p2_x + 1.2sphere_radius, 0.0, 0.0, -0.1, 0.0, 0.0]
    enter_problem = ODEProblem(constant_velocity!, enter_state, (0.0, 10.0), aux)
    enter_solution = solve(
        enter_problem,
        Vern7();
        callback=cb_enter(enter_events),
        abstol=1e-12,
        reltol=1e-12,
    )
    @test length(enter_events) == 1
    @test enter_events[1].code == :enter
    @test SEMDynamics.Dynamics.state_from_r1_r2(enter_events[1].state, μ)[2] ≈ sphere_radius atol=1e-10
    @test enter_solution.t[end] ≈ enter_events[1].time

    escape_events = dynamic_events()
    escape_state = [p2_x + 0.8sphere_radius, 0.0, 0.0, 0.1, 0.0, 0.0]
    escape_problem = ODEProblem(constant_velocity!, escape_state, (0.0, 10.0), aux)
    escape_solution = solve(
        escape_problem,
        Vern7();
        callback=cb_escape(escape_events),
        abstol=1e-12,
        reltol=1e-12,
    )
    @test length(escape_events) == 1
    @test escape_events[1].code == :escape
    @test SEMDynamics.Dynamics.state_from_r1_r2(escape_events[1].state, μ)[2] ≈ sphere_radius atol=1e-10
    @test escape_solution.t[end] ≈ escape_events[1].time
    @test_throws ArgumentError cb_enter(dynamic_events(); scale=0.0)
    @test_throws ArgumentError cb_escape(dynamic_events(); scale=-1.0)

    function run_apse_test(primary)
        primary_x = primary === :p1 ? -μ : 1 - μ
        center = primary_x + 0.2
        amplitude = 0.05
        initial_phase = 0.1

        function oscillator!(du, u, p, t)
            du[1] = u[3]
            du[2] = u[4]
            du[3] = -(u[1] - center)
            du[4] = -u[2]
            return nothing
        end

        state = [
            center + amplitude * cos(initial_phase),
            0.0,
            -amplitude * sin(initial_phase),
            0.0,
        ]
        events = dynamic_events()
        callback = if primary === :p1
            cb_apse_p1(
                events;
                perigee_scale=(0.14, 0.16),
                apogee_scale=(0.24, 0.26),
            )
        else
            cb_apse_p2(
                events;
                perilune_scale=(0.14, 0.16),
                apolune_scale=(0.24, 0.26),
            )
        end
        problem = ODEProblem(oscillator!, state, (initial_phase, 2pi + 0.1), aux)
        solve(problem, Vern7(); callback, abstol=1e-12, reltol=1e-12)
        return events
    end

    p1_events = run_apse_test(:p1)
    @test getproperty.(p1_events, :code) == [:perigee, :apogee]
    @test getproperty.(p1_events, :time) ≈ [pi, 2pi] atol=1e-10

    p2_events = run_apse_test(:p2)
    @test getproperty.(p2_events, :code) == [:perilune, :apolune]
    @test getproperty.(p2_events, :time) ≈ [pi, 2pi] atol=1e-10
end
