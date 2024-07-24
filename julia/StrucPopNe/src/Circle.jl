# Return a transition rate matrix describing two lineages within a
# circular stepping-stone model of population structure. The
# metapopulation consists of d demes. In state i, the two lineages are
# separated by i-1 steps. So state 1 means they are in the same deme,
# state 2 means they're separated by 1 steps, i.e. they are in
# adjacent demes. Time is measured in units of 2N generations. M = 4Nm
# is the instantaneous migration rate for a pairs of genes per unit of
# time. 
function rate_mat(::Val{:circle}, d::Integer, M::Real, notused::Real)
    maxdiff = d ÷ 2
    dim = maxdiff+1
    Q = zeros(Float64, dim, dim)
    Q[1,1] = -M - 1
    Q[1,2] = M
    for i in 2:maxdiff
        Q[i, i-1] = M/2
        Q[i, i+1] = M/2
        Q[i, i]   = -M
    end
    if iseven(d)
        Q[dim, dim-1] = M
        Q[dim, dim] = -M
    else
        Q[dim, dim-1] = M/2
        Q[dim, dim] = -M/2
    end
    return Q
end

# Calculate Fst for the circular stepping stone model of population
# structure. Approximation for small mutation rate. First presented on
# p 582 of Wilkinson-Herbots, Hilde M. (1998). Genealogy and
# Subpopulation Differentiation under Various Models of Population
# Structure Jol. of Math. Biol. 37(6): 535-585. doi:
# 10.1007/s002850050140. 
# - d: number of demes
# - M: migration rate, twice the number of immigrant genes per deme
#      per generation
function fst(::Val{:circle}, d::Integer, M::Real, notused::Real)
    h = 6*M/(d-(1/d))
    return 1/(1 + h)
end

# Wilkinson-Herbots's formula for global Fst in a circular stepping
# stone model. See Eqn. 34 of Wilkinson-Herbots, Hilde M. (1998).
# Genealogy and Subpopulation Differentiation under Various Models of
# Population Structure Jol. of Math. Biol. 37(6): 535-585. doi:
# 10.1007/s002850050140.
function fst_herbots(::Val{:circle}, d::Integer, M::Real, θ::Real)
    a = 2*π/d
    h = 0.0
    for k in 1:(d-1)
        h += 1/(θ + M*(1 - cos(a*k)))
    end
    return h/(d + h)
end

