using Revise
using BioprocessModels
using ModelingToolkit
using Latexify
using Symbolics
using Test

@named R = tank(nrInFlows=1, nrOutFlows=1)

# Define the components
@named X = component()
@named S = component(catalyst=X)
# @named P = component(catalyst=X)
# @named CO2 = component(catalyst=X)
# @named O2 = component(catalyst=X)
components = [X, S]#, P, CO2, O2]

# Define the reactions and couple to components
@named mnd = monod()
@named hld = haldane()
@named qiS = kinetic(mnd)
# @named mndhld = monod_haldane()
# @named qiP = kinetic(mndhld)
# @named mndhld2 = compose_kinetics(kin1=mnd, kin2=hld)
# @named qiP2 = kinetic(mndhld2)
reactions = [qiS]

# Define the reaction system
@named react = reaction_system(components, reactions)

# Compose 


full_equations(react)
parameters(react)
equations(qiP)
latexify(parameters(qiP2))
@named components = compose(R, X, S, P, CO2, O2)


@testset "BioprocessModels.jl" begin
    # Write your tests here.
end
