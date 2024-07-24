# Return a 2x2 transition rate matrix describing two lineages. In
# state 1, the 2 lineages are in a single deme, in state 2, they're in
# separate demes, in the absorbing state, 0, they have
# coalesced. Transition rates are 2->1: M/(d-1), 1->2: M, 1->0:
# 1. Here time is measured in units of 2N generations, and M = 4Nm,
# where m is the migration rate per generation.
function rate_mat(::Val{:island}, d::Integer, M::Real, notused::Real)
    Q = Matrix{Float64}(undef, 2, 2)
    Q[1,1] = -M - 1
    Q[1,2] = M
    Q[2,1] = M/(d-1)
    Q[2,2] = -M/(d-1)
    return Q
end

# Like rate_mat, but in discrete time.
function discrete_mat(::Val{:island}, d::Integer, m::Real, N::Real)
    p0 = (1-m)^2 # Pr neither gene migrates
    p1 = 2*m*(1-m) # Pr one migrates
    p2 = m^2 # Pr both migrate

    A = Matrix{Float64}(undef, 2, 2)
    A[1,1] = (p0 + p2/(d-1))*(1-1/(2N))
    A[1,2] = (p1 + p2*(d-2)/(d-1))*(1-1/(2N))
    A[2,1] = (p1+p2)/(d-1)
    A[2,2] = p0 + (p1+p2)*(d-2)/(d-1)
    return A
end

# Calculate Fst for the island model of population structure.
# - d: number of demes
# - M: equals 4Nm twice the number of immigrant genes per generation.
function fst(::Val{:island}, d::Integer, M::Real, ::Real)
    return 1/(1 + M*d^2/(d-1)^2)
end
