# Plotting, ex. 1
using Random: Xoshiro
rng = Xoshiro(1234)
xp1 = 0.2randn(rng, 10_000) .+ 0.5
p2 = 0.5randn(rng, 10_000)
edges = -2:0.1:2
hist(p1; bins = edges, color = ("black", 0.75), label = "p1")
hist!(p2; bins = edges, color = ("red", 0.75), label = "p2")
ylims!(0, nothing)
axislegend()
current_figure()