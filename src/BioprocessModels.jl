module BioprocessModels
export reaction_system, compose_kinetics, flow, gasflow, tank, component, kinetic, monod, blackman, teissier, moser, haldane, monod_haldane
using ModelingToolkit
using Symbolics

# Defining the independent variable time t and the differential operator D
@parameters t
D = Differential(t)

# Now define the components / subsystems of a bioprocess model
"""Defines a simple liquid flow rate from a reservoir into the tank or 
from a tank into a reservoir."""
function flow(; name)
    @parameters (
        ρ=1000, [description="density [g/L]"]
    )
    @variables (
        F(t), [description="flow rate [L/h]"],
        F_m(t), [description="flow rate [g/h]"]
        )
    eqs = [
        F_m ~ ρ * F,
    ]
    ODESystem(eqs, t, [F_m,F], [ρ], name = name)
end
# @named fl = flow()


function gasflow(; name)
    @constants (
        V_nM = 22.414, [description="molar volume [L/mol]"]
    )
    @variables (
        F(t), [description="gasflow rate [L/h]"],
        F_n(t), [description="gasflow rate [mol/h]"]
        )
    eqs = [
        F_n ~ F / V_nM,
    ]
    ODESystem(eqs, t, name = name)
end
# @named fng = gasflow()
"""
tank(; name, nrInFlows, nrOutFlows)

Defines a tank model with a specified name, number of inflows and number of outflows.

# Arguments
- `name::String`: Name of the tank model.
- `nrInFlows::Int`: Number of inflows to the tank.
- `nrOutFlows::Int`: Number of outflows from the tank.

# Returns
- `ODESystem`: An ODESystem object representing the tank model.
"""
function tank(; name, nrInFlows=2, nrOutFlows=1, nrGasFlows=1)
    @parameters (
        ρ=1000, [description="density [g/L]"]
    )
    @variables (
        V_L(t), [description="liquid volume [L]"],
        )
    inFlows = [flow(name=Symbol("inflow$i")) for i in 1:nrInFlows]
    outFlows = [flow(name=Symbol("outflow$i")) for i in 1:nrOutFlows]
    gasFlows = [gasflow(name=Symbol("gasflow$i")) for i in 1:nrGasFlows]
    eqs = [
        D(V_L) ~ sum([inFlow.F for inFlow in inFlows]) - sum([outFlow.F for outFlow in outFlows]),
        ]
    compose(ODESystem(eqs, t, [V_L], [ρ], name = name), inFlows..., outFlows...)
end
# @named F10 = tank()


function component(; name, catalyst=false)
    @parameters (
        M, [description="molar mass [g/mol]"]
    )
    @variables (
        m(t), [description="mass [g]"],
        n(t), [description="amount [mol]"],
        r(t), [description="reaction rate [g/h]"],
        r_n(t), [description="reaction rate [mol/h]"],
        q(t), [description="specific reaction rate [g/g/h]"],
        )
    eqs = [
        isa(catalyst,ModelingToolkit.AbstractODESystem) ? r ~ q * catalyst.m : r ~ q * m,
        r_n ~ r / M,
        n ~ m / M,
    ]
    ODESystem(eqs, t, [m,n,r,r_n,q], [M], name = name)
end

function reaction_system(componentsArray, reactionsArray; name)
    qi = [reaction.r for reaction in reactionsArray]
    q = [component.q for component in componentsArray]
    @parameters Y
    for (index_c, component) in enumerate(componentsArray)
        if index_c > 1
            Y = hcat(Y, [eval(Meta.parse("@parameters Y_" * String(component.name) * "_" * String(reaction.name)))[1] for reaction in reactionsArray])
        else
            Y = [eval(Meta.parse("@parameters Y_" * String(component.name) * "_" * String(reaction.name)))[1] for reaction in reactionsArray]
        end
    end
    Y = Y'
    @variables r(t), [description="Reaction rate"]
    #@parameters Y[1:length(componentsArray), 1:length(reactionsArray)], [description="Yields [g/g]"]
    compose(ODESystem(q .~ Y * qi, t, name = name), componentsArray..., reactionsArray...)
end

function compose_kinetics(; name, kin1::ModelingToolkit.AbstractODESystem, kin2::ModelingToolkit.AbstractODESystem)
    @variables r̃(t), [description="Reaction rate"]
    compose(ODESystem([r̃ ~ kin1.r̃ * kin2.r̃], name = name), kin1, kin2)
end
function monod_haldane(; name)
    @variables r̃(t), [description="Reaction rate"]
    @named mnd = monod()
    @named hld = haldane()
    compose(ODESystem([r̃ ~ mnd.r̃ * hld.r̃], name = name), mnd, hld)
end
"denormalize the kinetic expressions"
function kinetic(kin::ModelingToolkit.AbstractODESystem; name, r_max=1)
    @parameters r_max=r_max, [description="Maximal reaction rate"]
    @variables r(t), [description="Reaction rate"]
    @unpack r̃ = kin
    extend(ODESystem([r ~ r_max * r̃], name = name), kin)
end

function monod(; name, K=0.05)
    @parameters K=K, [description="Monod constant"]
    @variables (
        r̃(t), [description="Reaction rate"],
        c(t), [description="Substrate concentration"]
        )
    ODESystem([r̃ ~ c / (K + c)], t, [r̃,c], [K], name = name)
end
# @named black = blackman(K=0.1)
# @named qiS = kinetic(black, r_max = 1.2)

function blackman(; name, K=0.05)
    @parameters K=K, [description="Blackman constant"]
    @variables (
        r̃(t), [description="Reaction rate"],
        c(t), [description="Substrate concentration"]
        )
    ODESystem([r̃ ~ min(1, K * c)], t, [r̃,c], [K], name = name)
end

function teissier(; name, K=0.05)
    @parameters K=K, [description="Teissier constant"]
    @variables (
        r̃(t), [description="Reaction rate"],
        c(t), [description="Substrate concentration"]
        )
    ODESystem([r̃ ~ 1 - exp(-K * c)], t, [r̃,c], [K], name = name)
end

function moser(; name,K=0.05,N=2)
    @parameters (
        K=K, [description="Moser constant"],
        N, [description="Moser order"]
    )
    @variables (
        r̃(t), [description="Reaction rate"],
        c(t), [description="Substrate concentration"]
        )
    ODESystem([r̃ ~ c^N / (K^N + c^N)], t, [r̃,c], [K,N], name = name)
end

"haldane kinetics"
function haldane(; name, K=0.05, K_i=0.05)
    @parameters (
        K=K, [description="Haldane constant"],
        K_i=K_i, [description="Haldane inhibition constant"]
    )
    @variables (
        r̃(t), [description="Reaction rate"],
        c(t), [description="Substrate concentration"]
        )
    ODESystem([r̃ ~ c / (K + c * (1 + c / K_i))], t, [r̃,c], [K,K_i], name = name)
end

# macro species(name)
#     component(name = name)
# end

#@species T1

####### TESTS #######
# First build the components vector 

# X = component(name="X")
# component.(name=namevector, catalyst = X)
# @named S = component(catalyst = X)
# Y = ...
# components * 

end
