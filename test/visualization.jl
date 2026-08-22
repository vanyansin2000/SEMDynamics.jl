@testset "CR3BP scene setup" begin
    figure2d = Figure()
    axis2d = Axis(figure2d[1, 1])
    @test cr3bp_set_scn!(axis2d) === axis2d
    @test_throws ArgumentError cr3bp_set_scn!(axis2d; lims=(0.0, 1.0, 0.0))

    figure3d = Figure()
    axis3d = Axis3(figure3d[1, 1])
    @test cr3bp_set_scn!(axis3d) === axis3d
    @test cr3bp_set_scn!(axis3d; lims=(0.4, 1.4, -0.4, 0.4, -0.2, 0.2)) === axis3d
    @test_throws ArgumentError cr3bp_set_scn!(axis3d; lims=(0.0, 1.0, 0.0, 1.0))
end
