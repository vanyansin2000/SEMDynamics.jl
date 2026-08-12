using CairoMakie
using SEMDynamics

# 2:1 平面 DRO；初值、周期和 Jacobi 常数分别在 `.x0`、`.P` 和 `.C` 中。
dro_2_1 = generate_DRO(P=pi)

dro_3_1 = generate_DRO(P=pi * 2 / 3)

dro_4_1 = generate_DRO(P=pi * 2 / 4)

figure = Figure(size = (400 , 300))
axis = Axis(figure[1, 1], aspect=AxisAspect(1))
cr3bp_set_scn!(axis)
orb1 = plotPO!(axis, dro_2_1; linewidth=2)
orb2 = plotPO!(axis, dro_3_1; linewidth=2)
orb3 = plotPO!(axis, dro_4_1; linewidth=2)
axislegend(axis, [orb1 , orb2 , orb3] , ["2:1 DRO", "3:1 DRO", "4:1 DRO"] , position=:rt)
display(figure)
