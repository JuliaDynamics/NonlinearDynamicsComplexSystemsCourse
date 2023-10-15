# %% basins of attraction of multistable predator prey
using DynamicalSystems

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

# %%
# Decide how to plot attractors: go from attractor to real number
using Statistics: mean
a2r = A -> mean(x[1] for x in A)

plot_basins_attractors_curves(fractions_curves, attractors_info, a2r, prange)

# Now analyze to find when the limit cycle becomes large enough
# to be detected by our grid resolution
id = 2 # id of the periodic attractor (from plot)
j = findfirst(atts -> length(atts[id]) > 1, attractors_info)
p_hopf = prange[j]

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

thomas = CoupledODEs(thomas_rule, ones(3), [0.16])
xg = yg = zg = range(-6.0, 6.0; length = 101)

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
mean_norm = mean(shock_norms)

scatter(columns(A)...; color = shock_norms)
