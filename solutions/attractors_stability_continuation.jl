# %% basins of attraction of multistable predator prey
function predator_prey_rule(u, p, t)
    r, c, μ, ν, α, β, χ, δ = p
    N, P = u
    common = α*N*P/(β+N)
    dN = r*N*(1 - (c/r)*N)*((N-μ)/(N+ν)) - common
    dP = χ*common - δ*P
    return SVector(dN, dP)
end

u0 = SVector(8.0, 0.01)
r1, r2 = 1.0, 2.0
# r, c, μ, ν, α, β, χ, δ = p
p = [r1, 0.19, 0.03, 0.003, 800, 1.5, 0.004, 2.2]

using OrdinaryDiffEq: Tsit5
diffeq = (alg = Tsit5(), adaptive = false, dt = 0.05)
ds = CoupledODEs(predator_prey_rule, u0, p; diffeq)

density = 201
xg = range(0, 20; length = density)
yg = range(0, 0.03; length = density)
grid = (xg, yg)

colors = COLORS = [
    "#7143E0",
    "#191E44",
    "#0A9A84",
]

fig = Figure(resolution = (800, 400))
ax1 = Axis(fig[1,1]; title = "r = $(r1)", xlabel = "x", ylabel = "y")
ax2 = Axis(fig[1,2]; title = "r = $(r2)", xlabel = "x", yticklabelsvisible = false)

for (j, r) in enumerate((r1, r2))
    set_parameter!(ds, 1, r)
    mapper = AttractorsViaRecurrences(
        ds, grid; mx_chk_fnd_att = 2000, mx_chk_loc_att = 4000
    )
    ax = (ax1, ax2)[j]
    basins, attractors = basins_of_attraction(mapper, grid)
    heatmap!(ax, xg, yg, basins; colormap = colors)
    for k ∈ keys(attractors)
        x, y = columns(attractors[k])
        scatter!(ax, vec(attractors[k]);
            color = COLORS[k], markersize = 20,
            strokewidth = 3, strokecolor = :white
        )
    end
end

fig

# %% continuation of multistable predator prey
mapper = AttractorsViaRecurrences(
    ds, grid; mx_chk_fnd_att = 2000, mx_chk_loc_att = 4000
)

sampler, = statespace_sampler(
    min_bounds = minimum.(grid), max_bounds = maximum.(grid),
)

rsc = RecurrencesSeededContinuation(mapper)

pidx = 1 # index of parameter to change (here Kd)
prange = range(1, 2; length = 101) # parameter range to scan

fractions_curves, attractors_info = continuation(
    rsc, prange, pidx, sampler;
    show_progress = true, samples_per_parameter = 100
)

include(
    joinpath(pathof(Attractors), "..", "..", "docs", "basins_plotting.jl")
)

# Decide how to plot attractors: go from attractor to real number
using Statistics: mean
attractor_to_real = A -> mean((x[1] + x[2])/2 for x in A)

basins_attractors_curves_plot(fractions_curves, attractors_info, attractor_to_real, prange)
