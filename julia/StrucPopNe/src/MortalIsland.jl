# Mortal-island model: the finite island model with random extinctions.

# Return a 2x2 transition rate matrix describing two lineages. In
# state 1, the 2 lineages are in a single deme, in state 2, they're in
# separate demes, in the absorbing state, 0, they have
# coalesced. Transition rates are 2->1: (M+X)/(d-1), 1->2: M, 1->0:
# 1. Here time is measured in units of 2N generations, and M = 4Nm,
# where m is the migration rate per generation, and X = 4Nx, where x
# is the extinction rate per generation.
function rate_mat(::Val{:xisland}, d::Integer, M::Real, X::Real)
    Q = Matrix{Float64}(undef, 2, 2)
    Q[1,1] = -M - 1
    Q[1,2] = M
    Q[2,1] = (M+X)/(d-1)
    Q[2,2] = -M/(d-1)
    return Q
end

# Calculate Fst for the island model of population structure.
# - d: number of demes
# - M/2: rate of migration per 2N generations
# -X/2: rate of extinction per 2N generations
function fst(::Val{:xisland}, d::Integer, M::Real, X::Real)
    return 1/((M*d + X)*d/(d-1)^2 + 1)
end

# Vector of probabilities of states for a pair of genes chosen at
# random from the entire population.
function pr_t(::Union{Val{:island}, Val{:xisland}}, d::Integer)
    return [1/d, 1-1/d]
end
