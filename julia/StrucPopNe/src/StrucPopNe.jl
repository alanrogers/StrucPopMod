module StrucPopNe

using Gadfly, Printf, LinearAlgebra
export rate_mat, discrete_mat, rhaz, mkplot, fst, fst_herbots, pr_t
export ratio_of_change, get_d

include("Island.jl")
include("MortalIsland.jl")
include("Circle.jl")
include("MortalCircle.jl")

# Relative hazard. Full hazard of a coalescent event at time t
# is rhaz(Q, t)/(2N). Ne is N/rhaz(A,t).
function rhaz(Q::Matrix{Float64}, t::Real)
    v = exp(Q*t)[1,:]
    return v[1]/sum(v)
end

# Find d from heterozygosity and other parameters. Solves
# h = 4*N*u*(d+a)/(1+a) for d.
# - h: heterozygosity
# - u: mutation rate per site per generation
# - N: deme size
# - a: ratio of extinction rate to migration rate
function get_d(h, u, N, a)
    r = h/(4*N*u)
    return r*(1+a) - a
end

# tmax is maximum time in thousands of years.
# Xover M is X/M; extinction relative to migration
function mkplot(modsym::Symbol, d::Integer, M::Vector{T}, XoverM::Real,
                N::Real, tmax::Real;
                showNd=true,
                title::String = "") where {T <: Real}
    model = Val{modsym}()

    # convert time in units of thousands of years into units
    # of 2N generations, assuming generations of 29 years.
    scale_time(t) = (1000*t/29)/(2N)
    
    Q = [rate_mat(model, d, m, m*XoverM) for m in M]
    Nₑ = [t->N/rhaz(x,scale_time(t)) for x in Q] # vector of functions of t

    plt = plot(Nₑ, 0, tmax,
               Guide.colorkey(title="M", labels=string.(M)),
               Guide.xlabel("Years Ago / 1000"),
               Guide.ylabel("Nₑ"),
               Theme(minor_label_font_size=12pt, major_label_font_size=14pt,
                     key_label_font_size=12pt,
                     line_width=1pt)
               )
    if showNd
        ndline = layer(
            yintercept=[N*d],
            Geom.hline(color="gray")
        )
        plt = push!(plt, ndline)
    end
    if length(title) > 0
        plt = push!(plt, Guide.title(title))
    end
    return plt
end

# Ratio of asymptotic Ne to initial Ne. I.e., Ne(infinity)/Ne(0).
function ratio_of_change(model, d::Integer, M::Real, X::Real)
    M = rate_mat(model, d, M, X)
    v = inv(eigvecs(M))[end,:]
    return sum(v)/v[1]
end

end # module StrucPopNe
