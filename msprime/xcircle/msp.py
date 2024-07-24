#!/usr/bin/python3

import msprime as msp
import os, sys, textwrap
from random import expovariate, sample, randrange, random

def usage():
    print("Usage: ./msp.py [options]")
    print("  where options may include:")
    print("  -r or --run: run simulation. Default: run")
    print("               DemographyDebugger")
    sys.exit(1)

do_simulation = False
for arg in sys.argv[1:]:
    if arg == "-r" or arg == "--run":
        do_simulation = True
    else:
        usage()

# data
ohet = 2e-4 # observed heterozygosity in archaic genomes.
u_per_site = 1.4e-8 # mutation

# parameters
d = 11
N = 600 # deme size
M = 20 # rate of migration per pair of lineages per 2N generations
X = M # rate of extinction per pair of demes per 2N generations
m = M/(4*N) # migration rate per lineage per generation
x = X/(4*N) # extinction rate per deme per generation
max_t = 200000 # number of generations of simulated extinctions

nchromosomes = 100   # number of chromosomes
basepairs = 1e7     # number of nucleotides per chromosome
c = 1e-8 # recombination rate per base pair per generation

dem = msp.Demography()

# Populations
for k in range(d):
    dem.add_population(
        initial_size = N
    )

# Circular stepping-stone migration
for i in range(d):
    dem.set_migration_rate(i, (i+1)%d, m/2)
    dem.set_migration_rate(i, (i-1)%d, m/2)

# Extinction
global_extinction_rate = d*x

t = 0.0
n_extinctions = 0
while True:
    t += expovariate(global_extinction_rate)
    if t > max_t:
        break
    i = randrange(d)
    if random() < 0.5:
        dem.add_mass_migration(t, source=i, dest=(i+1)%d, proportion=1)
    else:
        dem.add_mass_migration(t, source=i, dest=(i-1)%d, proportion=1)
    n_extinctions += 1

print("n_extinctions =", n_extinctions, file=sys.stderr)

dem.sort_events()

if dem == None:
    raise ValueError("dem = None after sort_evenrts")
    
# One diploid sample from population 0
samples = [
    msp.SampleSet(num_samples=1, population=0, time=0, ploidy=2)
]

if do_simulation:
    heterozygosity = 0
    # Seed for random number generator. Uses 4 bytes (32 bits) from
    # /dev/urandom.
    seed = int.from_bytes(os.urandom(4), "big")

    # Simulate gene genealogy. When sim_ancestry is called without
    # num_replicates, it returns a value of class
    # tskit.trees.TreeSequence. But when it's called with
    # num_replicates, it returns a "generator", which can be
    # used in a loop to deal with each replicate in turn.
    chr_generator = msp.sim_ancestry(samples = samples,
                                     demography = dem,
                                     sequence_length = basepairs,
                                     random_seed = seed,
                                     recombination_rate = c,
                                     num_replicates = nchromosomes
    )
    
    # Put mutations onto the gene genealogy and output in psmcfa format.
    w = 100
    for i, chr in enumerate(chr_generator):
        #np.any(ts.nodes_time > demographic_model_time_limit)
        t = max(chr.nodes_time)
        if t > max_t:
            raise RuntimeError(f"t (={t}) > max_t (={max_t})")
        
        sim = msp.sim_mutations(chr, rate = u_per_site)
        L = int(sim.get_sequence_length() // w)
        outstr = ["T"] * L
        for v in sim.variants():
            gt = v.genotypes
            if gt[0] != gt[1]:
                heterozygosity += 1
                outstr[int(v.position // w)] = "K"
        print("> chr%d" % i)
        print("\n".join(textwrap.wrap("".join(outstr), width=79)))
        print()
    heterozygosity /= nchromosomes * basepairs
    print("heterozygosity: %0.5e" % heterozygosity, file=sys.stderr)
else:
    # Run demography debugger and quit
    print(dem.debug())
    print("Use \"./sim -r\" to run simulation")
