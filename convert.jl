# Convert tutorials of packages into Jupyter notebooks
import Pkg
Pkg.activate(@__DIR__)

import Literate
# Packages to convert tutorials
import DynamicalSystems, Attractors
# 3 tuple: package, name of tutorial file, name of output file
convert = [
    # (DynamicalSystems, "tutorial.jl", "dynamicalsystems_intro"),
    (Attractors, "tutorial.jl", "multiglobal_stability"),
]

for (modu, file, output) in convert
    inputfile = joinpath(pkgdir(modu), "docs", "src", file)
    Literate.notebook(
        inputfile, joinpath(@__DIR__, "notebooks");
        execute = false, credit = false, name = output,
    )
end
