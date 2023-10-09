
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
    # Estimation fails if `mx_chk_fnd_att` is too low, e.g., 10
    # this happens because near y=0 the trajectory slows down very much.
    # For this Δt value, this means that several integration steps are taken within
    # the _same_ cell of the grid. Hence, the algorithm counts this as recurrences.
    # 10 recurrences are very easy to accumulate in this way, and when this happens,
    # the algorithm switches to "new attractor found" mode and proceeds to (incorrectly)
    # identify some random cells near y = 0 as attractor cells.
    mx_chk_fnd_att = 100, Δt = 0.01,
)

basins, attractors = basins_of_attraction(mapper, basinsgrid)

fig = heatmap_basins_attractors(basinsgrid, basins, attractors)

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
