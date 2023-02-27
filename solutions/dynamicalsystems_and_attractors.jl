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
set_parameter!(ds, 1, 2.5)
for u in u0s
    X, t = trajectory(ds, 100.0, u)
    lines!(ax, t, X[:, 1])
end
fig