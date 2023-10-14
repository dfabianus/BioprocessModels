using Revise
using BioprocessModels
using ModelingToolkit
using Latexify
using Symbolics
using DifferentialEquations
using Plots
using Test

@named F1 = flow()
@named Fout = flow()
@named R = tank(F1)

# Define the components
@named X = component(R)
@named S = component(R,catalyst=X)
components = [X, S]

# Define the reactions and couple to components
@named mnd = monod(S)
@named qiS = kinetic(mnd)
reactions = [qiS]

# Define the reaction system
@named react = reaction_system(components, reactions)

@named feedsys = feed_system(components, [F1])

# Compose bioprocess system
@named bp = bioprocess(R, react, components)
bpsimp = structural_simplify(bp)

x0 = [
    R.V_L => 1.0,
    X.m => 0.1,
    S.m => 10.0,
]
ts = (0.0, 12.0)
p = [
    react.Y_X_qiS => 0.45,
    react.Y_S_qiS => -1.0,
    X.M => 26.5,
    S.M => 30.0,
]
prob = ODEProblem(bpsimp, x0, ts, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
plot(sol,idx=[X.m, S.m])

@testset "BioprocessModels.jl" begin
    # Write your tests here.
end
