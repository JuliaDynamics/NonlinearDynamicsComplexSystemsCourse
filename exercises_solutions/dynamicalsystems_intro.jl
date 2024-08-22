import Pkg; Pkg.activate(dirname(@__DIR__))
using DynamicalSystems, CairoMakie

# Exercise 1
function roessler_rule(u, p, t)
    x, y, z = u
    a, b, c = p
    dx = -y - z
    dy = x + a*y
    dz = b + z*(x - c)
    return SVector(dx, dy, dz)
end

u0 = [1, -2, 0.1]
p0 = [0.2, 0.2, 5.0]
roessler = CoupledODEs(roessler_rule, u0, p0)

X, t = trajectory(roessler, 1000.0; Ttr = 100, Dt = 0.01)
lines(columns(X)...; axis = (type = Axis3,))

# %% Exercise 2
as = range(0.1, 0.3; length = 101)
cs = range(1, 10; length = 101)

function lyapunov_from_params(a, c)
    set_parameter!(roessler, 1, a)
    set_parameter!(roessler, 3, c)
    reinit!(roessler)
    return lyapunov(roessler, 1000.0; Ttr = 100.0)
end

# Obtian the exponents
@time ls = broadcast(lyapunov_from_params, as, cs')
# Find parameters with positive
positive_idxs = findall(λ -> abs(λ) > 1e-3, ls)
params = [(a = as[c[1]], c = cs[c[2]]) for c in positive_idxs]
# Fraction
chaotic_fraction = length(params)/length(ls)

heatmap(as, cs, ls)

# %% Exercise 3
@inbounds function ikedamap_rule(u, p, n)
    a,b,c,d  = p
    t = c - d/(1 + u[1]^2 + u[2]^2)
    dx = a + b*(u[1]*cos(t) - u[2]*sin(t))
    dy = b*( u[1]*sin(t) + u[2]*cos(t) )
    return SVector(dx, dy)
end

u0 = [1.0, 1.0]
p1 = (a=1.0, b=1.0, c=0.4, d =6.0)
p2 = [6, 0.9, 3.1, 6]
ikedamap = DeterministicIteratedMap(ikedamap_rule, u0, p2)

X, t = trajectory(ikedamap, 10_000; Ttr = 100)

fig, = scatter(columns(X)...)

@show grassberger_proccacia_dim(X)

X1 = StateSpaceSet([p for p in X if p[1] < 2.5])
X2 = StateSpaceSet([p for p in X if p[1] > 2.5])

@show grassberger_proccacia_dim(X1)
@show grassberger_proccacia_dim(X2)


# %% Poincare tristable

function thomas_rule(u, p, t)
    x,y,z = u
    b = p[1]
    xdot = sin(y) - b*x
    ydot = sin(z) - b*y
    zdot = sin(x) - b*z
    return SVector(xdot, ydot, zdot)
end

thomas = CoupledODEs(thomas_rule, 3rand(3), [0.1665])

pmap = PoincareMap(thomas, (3, 0.0))

X, t = trajectory(pmap, 1000)
X

# Plot many initial conditions
    fig = Figure(size = (1000, 400))
for (i, b) in enumerate((0.16, 0.17, 0.18))
    set_parameter!(pmap, 1, b)
    ax = Axis(fig[1,i]; limits = ((-6, 6), (-6, 6)))
    for j in 1:20
        reinit!(pmap, 3rand(3))
        X, t = trajectory(pmap, 10; Ttr = 100)
        scatter!(ax, X[:, 1], X[:, 2]; color = Cycled(j), markersize = rand()*15 + 5)
    end
end

fig


# %%

reinit!(pmap, 3rand(3))
current_state(pmap)
X, t = trajectory(pmap, 10; Ttr = 100)
X[end]