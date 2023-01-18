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

# Plotting, ex. 2
fig = Figure(resolution = (600, 600))

colors = [
    "#7143E0",
    "#191E44",
    "#0A9A84",
    "#C0A12B",
    "#701B80",
    "#2E6137",
]

Box(fig[1, 1:2], color = colors[1], strokewidth = 0)
Box(fig[1:2, 3], color = colors[2], strokewidth = 0)
Box(fig[3, 2:3], color = colors[3], strokewidth = 0)
Box(fig[2:3, 1], color = colors[4], strokewidth = 0)
Box(fig[2, 2], color = colors[5], strokewidth = 0)

ax = Axis(fig[:, :])
n = 4
t = 0:0.01:2n*π
x = (2n*π .- t) .* cos.(t)
y = (2n*π .- t) .* sin.(t)

lines!(ax, x, y; color = colors[6], linewidth = 5)
hidedecorations!(ax)
hidespines!(ax)
fig
