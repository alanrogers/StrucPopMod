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

# Parameters. We begin with M = 4Nm, where N is deme size and m is the
# rate of migration per lineage per generation in backwards continuous time. 
M = 10
d = 6  # number of demes
N = 600 # size of each deme
m = M/(4*N)

nchromosomes = 10   # number of chromosomes
basepairs = 1e7     # number of nucleotides per chromosome
u_per_site = 1.4e-8 # mutation
c = 1e-8 # recombination rate per base pair per generation

# Circular stepping-stone model. First argument is a list of d deme
# sizes, each of which equals N. Second argument is the rate of
# migration between each pair of neighboring demes. Third argument
# equals False, which means that the two ends of the sequence of demes
# are joined to form a circle. In other words, there are no
# boundaries.
dem = msp.Demography.stepping_stone_model(d*[N], m/2, boundaries=False)

if dem == None:
    raise ValueError("dem = None after Demography")
    
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
    print("Use \"python msp.py -r\" to run simulation")
