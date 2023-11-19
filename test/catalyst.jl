using Catalyst
rn = @reaction_network begin
    k1, 6C --> X
end

reactions(rn)
Graph(rn)
equations(rn)
odesys = convert(ODESystem, rn)
equations(odesys)

