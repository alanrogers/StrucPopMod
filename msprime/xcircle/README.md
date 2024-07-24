# Simulations under circular stepping-stone model with extinction

This model assumes *d* demes of size *N*. In backwards time, each pair
of demes exchanges migrants at rate *M* per unit of *2N*
generations. In addition, demes may also go extinct. In forward time,
extinction means that a random deme disappears and is immediately
replaced by immigrants from a single other deme, which is chosen at
random. In backwards time, extinction events look like mass migration,
with all lineages in the affected deme moving simultaneously to a
different (randomly-chosen) deme.


    python3 msp.py -r 1>sim1.psmcfa 2>heterozygosity.txt
    python3 msp.py -r 1>sim2.psmcfa 2>>heterozygosity.txt
    python3 msp.py -r 1>sim3.psmcfa 2>>heterozygosity.txt

    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim1.psmc sim1.psmcfa &
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim2.psmc sim2.psmcfa &
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim3.psmc sim3.psmcfa &

    psmcdata -u 1.4e-8 -g 29 sim1.psmc > sim1.txt
    psmcdata -u 1.4e-8 -g 29 sim2.psmc > sim2.txt
    psmcdata -u 1.4e-8 -g 29 sim3.psmc > sim3.txt

Then within julia:

    include("simplot.jl")
	savefig("simxcircle.png")



