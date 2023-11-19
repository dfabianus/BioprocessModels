using Catalyst

rn = @reaction_network begin
    yXS * qSmax * S / (kS + S), X --> X
    qSmax, S --> 0
end

reactions(rn)
Graph(rn)
equations(rn)
odesys = convert(ODESystem, rn)
equations(odesys)

