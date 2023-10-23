module BioprocessModels
export feed_system, bioprocess, reaction_system, compose_kinetics, flow, gasflow, tank, component, kinetic, monod, blackman, teissier, moser, haldane, monod_haldane
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
    ρ = DelayParentScope(ρ)
    F = DelayParentScope(F)
    F_m = DelayParentScope(F_m)
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
function tank(flows...; name)
    @parameters (
        ρ=1000, [description="density [g/L]"]
    )
    @variables (
        V_L(t), [description="liquid volume [L]"],
        )
    V_L = ParentScope(V_L)
    #inFlows = [flow(name=Symbol("inflow$i")) for i in 1:nrInFlows]
    #outFlows = [flow(name=Symbol("outflow$i")) for i in 1:nrOutFlows]
    #gasFlows = [gasflow(name=Symbol("gasflow$i")) for i in 1:nrGasFlows]
    eqs = [
        D(V_L) ~ sum([flow.F for flow in flows]),
        ]
    compose(ODESystem(eqs, t, [V_L], [ρ], name = name), flows...)
end
# @named F10 = tank()


function component(reactor; name, catalyst=false)
    @parameters (
        M, [description="molar mass [g/mol]"]
    )
    @variables (
        m(t), [description="mass [g]"],
        c(t), [description="concentration [g/L]"],
        Q_in(t), [description="incoming mass flow rate [g/h]"],
        n(t), [description="amount [mol]"],
        r(t), [description="reaction rate [g/h]"],
        r_n(t), [description="reaction rate [mol/h]"],
        q(t), [description="specific reaction rate [g/g/h]"],
        )
    V_L = ParentScope(reactor.V_L)
    eqs = [
        isa(catalyst,ModelingToolkit.AbstractODESystem) ? r ~ q * ParentScope(catalyst.m) : r ~ q * m,
        r_n ~ r / M,
        n ~ m / M,
        c ~ m / V_L,
    ]
    ODESystem(eqs, t, [V_L,m,c,n,r,r_n,q,Q_in], [M], name = name)
end

function feed_system(componentsArray, incomingFeedsArray; name)
    @parameters c_F, [description="Feed concentration [g/L]"]
    @variables Q_in(t)
    for (index_c, component) in enumerate(componentsArray)
        if index_c > 1
            c_F = hcat(c_F, [eval(Meta.parse("@parameters c_" * String(component.name) * "_" * String(feed.name)))[1] for feed in incomingFeedsArray])
        else
            c_F = [eval(Meta.parse("@parameters c_" * String(component.name) * "_" * String(feed.name)))[1] for feed in incomingFeedsArray]
        end
    end
    c_F = ParentScope.(c_F)
    c_F = c_F'
    
    eqs = [
        [ParentScope(component.Q_in) for component in componentsArray] .~ c_F * [ParentScope(feed.F) for feed in incomingFeedsArray]
    ]
    ODESystem(reduce(vcat,eqs), t, name = name)
end


function bioprocess(reactor, reactionSystem, componentsArray, feedSystem; name)
    @unpack F1 = reactor
    eqs = [
        [D(component.m) for component in componentsArray] .~ [component.r for component in componentsArray] .+ [component.Q_in for component in componentsArray],
        F1.F ~ 0,
    ]
    compose(ODESystem(reduce(vcat,eqs), t, name = name), feedSystem, reactor, reactionSystem, componentsArray...)
end

function reaction_system(componentsArray, reactionsArray; name)
    qi = [reaction.r for reaction in reactionsArray]
    q = [component.q for component in componentsArray]
    q = ParentScope.(q)
    qi = ParentScope.(qi)
    @parameters Y
    for (index_c, component) in enumerate(componentsArray)
        if index_c > 1
            Y = hcat(Y, [eval(Meta.parse("@parameters Y_" * String(component.name) * "_" * String(reaction.name)))[1] for reaction in reactionsArray])
        else
            Y = [eval(Meta.parse("@parameters Y_" * String(component.name) * "_" * String(reaction.name)))[1] for reaction in reactionsArray]
        end
    end
    Y = ParentScope.(Y)
    Y = Y'
    @variables r(t), [description="Reaction rate"]
    #@parameters Y[1:length(componentsArray), 1:length(reactionsArray)], [description="Yields [g/g]"]
    compose(ODESystem(q .~ Y * qi, t, name = name), reactionsArray...)
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
    r = DelayParentScope(r)
    r_max = DelayParentScope(r_max)
    @unpack r̃ = kin
    r̃ = DelayParentScope(r̃)
    extend(ODESystem([r ~ r_max * r̃], name = name), kin)
end

function monod(substrate; name, K=0.05)
    @parameters K=K, [description="Monod constant"]
    @variables (
        r̃(t), [description="Reaction rate"],
        #c(t), [description="Substrate concentration"]
        )
    r̃ = DelayParentScope(r̃)
    c = ParentScope(ParentScope(substrate.c))
    K = DelayParentScope(K)
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

####### CONTROL #######
function qS_control(;name)
    @variables (
        F()
    )
end

end
