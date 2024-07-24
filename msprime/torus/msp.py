#!/usr/bin/python3

# Before running this: "python3 -m pip install msprime"
import msprime as msp
import os, sys, textwrap

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

# Regardless of the pattern of migration or the sizes of local
# populations, within-group heterozygosity should equal 4Ndu, where N
# is deme size, d is the number of demes, and u is mutation rate. With
# the data values above, this gives Nd = 3571. If d were 4, then
# N=893. However, you can't build a torus with only 4 demes. The smallest
# torus that isn't just the island model has a 4x4 grid: 16 demes in all.
# This requires N = 223.

# parameters
d = 144
dim = 12  # sqr root of number of demes. Must be at least 4.
N = 25 # deme size
m = 0.1

nchromosomes = 100   # number of chromosomes
basepairs = 2e6     # number of nucleotides per chromosome
c = 1e-8 # recombination rate per base pair per generation

dem = msp.Demography()

# Populations
for k in range(d):
    dem.add_population(
        initial_size = N
    )

# Migration on a torus. Makes use of facts that for positive n,
# -1 % n == n-1, and n % n == 0
for i in range(dim):
    for j in range(dim):
        src = i + dim*j # index of population (i,j)
        dst = i + dim*((j+1) % dim)
        dem.set_migration_rate(src, dst, m/4)
        dst = i + dim*((j-1) % dim)
        dem.set_migration_rate(src, dst, m/4)
        dst = ((i+1) % dim) + dim*j
        dem.set_migration_rate(src, dst, m/4)
        dst = ((i-1) % dim) + dim*j
        dem.set_migration_rate(src, dst, m/4)

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
