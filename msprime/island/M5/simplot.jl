using Plots, DataFrames, CSV, LaTeXStrings

# Read file and convert. Ne becomes an estimate of numbers of individuals,
# and t becomes time in years
function getdf(filename)
    df = CSV.read(filename, header=3, DataFrame)
    return df
end

popnames = ["sim1", "sim2", "sim3"]

# A vector of data frames
dfs = [getdf(pname*".txt") for pname in popnames];

for df in dfs
    df.t ./= 1000
end

plot(dfs[1].t, dfs[1].Ne, label=popnames[1],
     xlims=(0,50),
     #ylims=(0, 2.5), 
     thickness_scaling=1.5, legend=:bottomright)
for i in 2:length(popnames)
    plot!(dfs[i].t, dfs[i].Ne, label=popnames[i],
          xlims=(0,50),
          #ylims=(0, 2.5)
          )
end
xlabel!("Years Ago / 1000")
ylabel!("Nₑ")
title!("Island: M=5, N=600, d=6")


