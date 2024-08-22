import Pkg; Pkg.activate(dirname(@__DIR__))

# %% basins of attraction of multistable predator prey
using DynamicalSystems, CairoMakie

function predator_prey_rule(u, p, t)
    r, c, μ, ν, α, β, χ, δ = p
    N, P = u
    common = α*N*P/(β+N)
    dN = r*N*(1 - (c/r)*N)*((N-μ)/(N+ν)) - common
    dP = χ*common - δ*P
    return SVector(dN, dP)
end

u0 = SVector(8.0, 0.01)
# r, c, μ, ν, α, β, χ, δ = p
p = [2.0, 0.19, 0.03, 0.003, 800, 1.5, 0.004, 2.2]

using OrdinaryDiffEq: Rodas5P
diffeq = (alg = Rodas5P(), abstol = 1e-9, rtol = 1e-9)
ds = CoupledODEs(predator_prey_rule, u0, p; diffeq)

density = 101
xg = range(-0.1, 20; length = density)
yg = range(-0.0001, 0.03; length = density)
grid = (xg, yg)

basinsgrid = (xg[2:end], yg[2:end])

mapper = AttractorsViaRecurrences(
    ds, grid;
    # Estimation fails if `consecutive_recurrences` is too low, e.g., 10
    # this happens because near y=0 the trajectory slows down very much.
    # For this Δt value, this means that several integration steps are taken within
    # the _same_ cell of the grid. Hence, the algorithm counts this as recurrences.
    # 10 recurrences are very easy to accumulate in this way, and when this happens,
    # the algorithm switches to "new attractor found" mode and proceeds to (incorrectly)
    # identify some random cells near y = 0 as attractor cells.
    consecutive_recurrences = 1000, Δt = 0.01,
)

basins, attractors = basins_of_attraction(mapper, basinsgrid)

fig = heatmap_basins_attractors(basinsgrid, basins, attractors)

fig

# %% continuation of multistable predator prey
sampler, = statespace_sampler(basinsgrid, 12345)

mapper = AttractorsViaRecurrences(
    ds, grid;
    consecutive_recurrences = 5000, Δt = 0.1,
    force_non_adaptive = true, Ttr = 100.0,
)

rafm = RecurrencesFindAndMatch(mapper)

pidx = 1 # index of parameter to change (here Kd)
prange = range(1, 2; length = 101) # parameter range to scan

fractions_curves, attractors_info = continuation(
    rafm, prange, pidx, sampler;
    show_progress = true, samples_per_parameter = 100
)

# Decide how to plot attractors: go from attractor to real number
using Statistics: mean
a2r = A -> mean(x[1] for x in A)

plot_basins_attractors_curves(fractions_curves, attractors_info, a2r, prange)

# Now analyze to find when the limit cycle becomes large enough
# to be detected by our grid resolution
id = 2 # id of the periodic attractor (from plot)
j = findfirst(atts -> length(atts[id]) > 1, attractors_info)
p_hopf = prange[j]


# %% Basin instability in predator pray (phase tipping)
using DynamicalSystems, CairoMakie

function predator_prey_rule(u, p, t)
    r, c, μ, ν, α, β, χ, δ = p
    N, P = u
    common = α*N*P/(β+N)
    dN = r*N*(1 - (c/r)*N)*((N-μ)/(N+ν)) - common
    dP = χ*common - δ*P
    return SVector(dN, dP)
end

u0 = SVector(8.0, 0.01)
# r, c, μ, ν, α, β, χ, δ = p
p = [2.0, 0.19, 0.03, 0.003, 800, 1.5, 0.004, 2.2]

using OrdinaryDiffEq: Rodas5P
diffeq = (alg = Rodas5P(), abstol = 1e-9, rtol = 1e-9)
ds = CoupledODEs(predator_prey_rule, u0, p; diffeq)

density = 101
xg = range(-0.1, 20; length = density)
yg = range(-0.0001, 0.03; length = density)
grid = (xg, yg)

basinsgrid = (xg[2:end], yg[2:end])

sampler, = statespace_sampler(basinsgrid, 12345)

r1 = 1.8
r2 = 2.5
set_parameter!(ds, 1, r1)

mapper = AttractorsViaRecurrences(
    ds, grid;
    consecutive_recurrences = 5000, Δt = 0.1,
    force_non_adaptive = true, Ttr = 100.0,
    sparse = false,
)

# basins at r1
# Plot basins before and after
fig = Figure(resolution = (800, 400))
basins = []
attractors = []
for i in 1:2
    r = (r1, r2)[i]
    set_parameter!(ds, 1, r)
    ax = Axis(fig[1,i]; title = "r = $(r)")
    basins1, attractors1 = basins_of_attraction(mapper, basinsgrid)
    heatmap_basins_attractors!(ax, basinsgrid, basins1, attractors1)
    push!(basins, copy(basins1))
    push!(attractors, copy(attractors1))
    Attractors.reset!(mapper)
end
fig

# Estimate basin instability
basins1, basins2 = basins
attractors1, attractors2 = attractors
# match:
rmap = match_statespacesets!(attractors2, attractors1)
replace!(basins2, rmap...)
# Estimate overlap basin that point in attractor goes to
# in fact, you do not need the basins. you can just evolve points
# on the new attractor... But anyways here we have the basins already.

grid_info = mapper.bsn_nfo.grid_nfo

basin_instabilities = Dict()

for (id, A) in attractors2
    # estimate overlap with basins1
    overlaps = zeros(Bool, length(A))
    for (i, p) in enumerate(A)
        cartesian_index = Attractors.basin_cell_index(p, grid_info)
        overlaps[i] = id == basins1[cartesian_index]
    end
    basin_instabilities[id] = overlaps
end

basin_instabilities[1]

# If this is true, there exists basin instability
any(isequal(false), basin_instabilities[1])
any(isequal(false), basin_instabilities[2])



# %% MFS of Thomas cyclical
using DynamicalSystems

function thomas_rule(u, p, t)
    x,y,z = u
    b = p[1]
    xdot = sin(y) - b*x
    ydot = sin(z) - b*y
    zdot = sin(x) - b*z
    return SVector(xdot, ydot, zdot)
end

thomas = CoupledODEs(thomas_rule, rand(3), [0.17])
xg = yg = zg = range(-6.0, 6.0; length = 100)

grid = (xg, yg, zg)

mapper = AttractorsViaRecurrences(thomas, grid)

id = mapper([1, 2, 3.0])
A = extract_attractors(mapper)[id]

mfs_algo = MFSBruteForce(1000, 1000, 0.99)
search_area = (-6.0, 6.0)
shocks = map(u0 -> minimal_fatal_shock(mapper, u0, search_area, mfs_algo), A)

using CairoMakie
using LinearAlgebra: norm
using Statistics: mean

shock_norms = norm.(shocks)
mean_shock = mean(shocks)
min_norm = minimum(shock_norms)
shock_norms ./= min_norm/10

fig = Figure(size = (800, 400))
for (i, a) in enumerate((0.7, 4.0))
    scatter(fig[1,i], columns(A)...; color = shock_norms,
        alpha = 0.5, axis = (type = Axis3, azimuth = a), markersize = 15,
    )
end
fig
