using GLMakie
using SEMDynamics

# Planar distant retrograde orbits.
dro_2_1 = generate_DRO(P=pi)
dro_3_1 = generate_DRO(P=2pi / 3)
dro_4_1 = generate_DRO(P=2pi / 4)

# Spatial L1 halo representatives from the northern and southern branches.
l1_northern = generate_halo(branch=:northern, lp=:L1 )
l1_southern = generate_halo(branch=:southern, lp=:L1)
l2_northern = generate_halo(branch=:northern, lp=:L2 )
l2_southern = generate_halo(branch=:southern, lp=:L2 )

# The northern L2 halo and the southern 9:2 NRHO representative.
northern_nrho_9_2 = generate_nrho_9_2(branch=:northern)
southern_nrho_9_2 = generate_nrho_9_2(branch=:southern)
figure = Figure(size=(1200, 420))

# DROs are planar, so a 2D axis shows their geometry most clearly.
dro_axis = Axis(
    figure[1, 1];
    title="Planar DROs",
    xlabel="x [LU]",
    ylabel="y [LU]",
    aspect=AxisAspect(1),
)
cr3bp_set_scn!(dro_axis)
dro_plots = [
    plotPO!(dro_axis, dro_2_1; color=:royalblue, linewidth=2),
    plotPO!(dro_axis, dro_3_1; color=:darkorange, linewidth=2),
    plotPO!(dro_axis, dro_4_1; color=:seagreen, linewidth=2),
]
axislegend(dro_axis, dro_plots, ["2:1 DRO", "3:1 DRO", "4:1 DRO"]; position=:rt)

# Northern and southern L1 halo branches.
l1_axis = Axis3(
    figure[1, 2];
    title="L1 halo branches",
    xlabel="x [LU]",
    ylabel="y [LU]",
    zlabel="z [LU]",
)
cr3bp_set_scn!(l1_axis; lims=(0.8, 1.2, -0.2, 0.2, -0.2, 0.2))
l1_plots = [
    plotPO!(l1_axis, l1_northern; color=:royalblue, linewidth=2),
    plotPO!(l1_axis, l1_southern; color=:darkorange, linewidth=2),
    plotPO!(l1_axis, l2_northern;  linewidth=2),
    plotPO!(l1_axis, l2_southern;  linewidth=2),
]
# axislegend(l1_axis, l1_plots, ["Northern L1 halo", "Southern L1 halo" , ""])

# The L2 panel contrasts the northern branch with the southern 9:2 NRHO.
l2_axis = Axis3(
    figure[1, 3];
    title="northern and southern 9:2 NRHO",
    xlabel="x [LU]",
    ylabel="y [LU]",
    zlabel="z [LU]",
)
cr3bp_set_scn!(l2_axis; lims=(0.8, 1.2, -0.2, 0.2, -0.2, 0.2))
l2_plots = [
    plotPO!(l2_axis, northern_nrho_9_2; color=:mediumpurple, linewidth=2),
    plotPO!(l2_axis, southern_nrho_9_2; color=:crimson, linewidth=2.5),
]
# axislegend(l2_axis, l2_plots, ["Northern L2 halo", "Southern 9:2 NRHO"])

display(figure)
