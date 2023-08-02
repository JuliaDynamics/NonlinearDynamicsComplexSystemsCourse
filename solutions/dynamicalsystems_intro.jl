using DynamicalSystems, CairoMakie
# Exercise 1

function predator_prey_rule(u, p, t)
    r, c, μ, ν, α, β, χ, δ = p
    N, P = u
    common = α*N*P/(β+N)
    dN = r*N*(1 - (c/r)*N)*((N-μ)/(N+ν)) - common
    dP = χ*common - δ*P
    return SVector(dN, dP)
end

u0 = SVector(8.0, 0.01)
r1, r2 = 1.7, 2.5
# r, c, μ, ν, α, β, χ, δ = p
p = [r1, 0.19, 0.03, 0.003, 800, 1.5, 0.004, 2.2]
ds = CoupledODEs(predator_prey_rule, u0, p)

u0s = [
    [5, 0.016],
    [5, 0.025],
    [15, 0.012],
]

fig = Figure()
ax = Axis(fig[1,1])
for u in u0s
    X, t = trajectory(ds, 100.0, u)
    lines!(ax, t, X[:, 1])
end

ax = Axis(fig[2,1])
set_parameter!(ds, 1, r2)
for u in u0s
    X, t = trajectory(ds, 100.0, u)
    lines!(ax, t, X[:, 1])
end
fig

# %% Exercise 2
using DynamicalSystems, GLMakie

@inbounds function ikedamap_rule(u, p, n)
    a,b,c,d  = p
    t = c - d/(1 + u[1]^2 + u[2]^2)
    dx = a + b*(u[1]*cos(t) - u[2]*sin(t))
    dy = b*( u[1]*sin(t) + u[2]*cos(t) )
    return SVector{2}(dx, dy)
end

u0 = [1.0, 1.0]
p1 = (a=1.0, b=1.0, c=0.4, d =6.0)
p2 = [6, 0.9, 3.1, 6]
ikedamap = DeterministicIteratedMap(ikedamap_rule, u0, p2)

X, t = trajectory(ikedamap, 10_000; Ttr = 100)

fig, = scatter(columns(X)...)
display(fig)

@show grassberger_proccacia_dim(X)

X1 = StateSpaceSet([p for p in X if p[1] < 2.5])
X2 = StateSpaceSet([p for p in X if p[1] > 2.5])

@show grassberger_proccacia_dim(X1)
@show grassberger_proccacia_dim(X2)
