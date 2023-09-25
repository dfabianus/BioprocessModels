module BioprocessModels

using ModelingToolkit

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

function gasflow(; name)
    @constants (
        V_nM = 22.4, [description="molar volume [L/mol]"]
    )
    @variables (
        F(t), [description="gasflow rate [L/h]"],
        F_n(t), [description="gasflow rate [mol/h]"]
        )
    eqs = [
        F_n ~ F / V_nM,
    ]
    ODESystem(eqs, t, [Q_L,F_L], [ρ], name = name)
end

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
    inFlows = [flow(name="inflow$i") for i in 1:nrInFlows]
    outFlows = [flow(name="outflow$i") for i in 1:nrOutFlows]
    gasFlows = @variables(F_G(t), :gasflow, 1:nrGasFlows)
    eqs = [
        D(V_L) ~ sum(F_L_in) - sum(F_L_out),

    ]
    compose(ODESystem(eqs, t, [V_L], [ρ], name = name), inFlows, outFlows)
end

function component(; name, catalyst::ModelingToolkit.AbstractODESystem=nothing)
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
        catalyst ? r ~ q * catalyst.m : r ~ q * m,
        r_n ~ r / M,
        n ~ m / M,
    ]
    ODESystem(eqs, t, [m,n,r,r_n,q], [M], name = name)
end

function monod(; name, K = 0.05, q_max = 1.1)
    
end

# X = component(name="X")
# component.(name=namevector, catalyst = X)
# @named S = component(catalyst = X)
# Y = ...
# components * 

end
