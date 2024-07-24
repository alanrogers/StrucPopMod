# Return a transition rate matrix describing two lineages within a
# circular stepping-stone model with extinction. The
# metapopulation consists of d demes. In state i, the two lineages are
# separated by i-1 steps. So state 1 means they are in the same deme,
# state 2 means they're separated by 1 steps, i.e. they are in
# adjacent demes.
# - d: number of demes
# - M: migration rate per pair of lineages per 2N generations
# - X: extinction rate per pair of deme per 2N generations.
function rate_mat(::Val{:xcircle}, d::Integer, M::Real, X::Real)
    maxdiff = d ÷ 2
    dim = maxdiff+1
    MpX = M+X
    Q = zeros(Float64, dim, dim)
    Q[1,1] = -M - 1
    Q[1,2] = M
    for i in 2:maxdiff
        Q[i, i-1] = MpX/2
        Q[i, i+1] = MpX/2
        Q[i, i]   = -MpX
    end
    if iseven(d)
        Q[dim, dim-1] = MpX
        Q[dim, dim] = -MpX
    else
        Q[dim, dim-1] = MpX/2
        Q[dim, dim] = -MpX/2
    end
    return Q
end

# Calculate Fst for the circular stepping stone model with extinction.
# - d: number of demes
# - M: migration rate per pair of lineages per 2N generations
# - X: extinction rate per pair of deme per 2N generations.
function fst(model::Val{:xcircle}, d::Integer, M::Real, X::Real)
    return 1.0/(6*(M*d + X)/(d^2-1) + 1.0)
end

# Vector of probabilities of states for a pair of genes chosen at
# random from the entire population.
function pr_t(::Union{Val{:circle}, Val{:xcircle}}, d::Integer)
    dim = (d ÷ 2) + 1
    v = Vector{Float64}(undef, dim)
    v[1] = 1/d
    if isodd(d)
        v[2:end] .= 2/d
    else
        v[2:end-1] .= 2/d
        v[end] = 1/d
    end
    return v
end


