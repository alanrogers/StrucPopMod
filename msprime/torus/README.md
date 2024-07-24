# Simulations under a model of isolation by distance on a torus

A torus is a geometrical object shaped like a doughnut. Start with a
square *L* by *L* lattice of demes, roll it up in one dimension to
form a cylinder with length *L* and circumference *L*. The join the
ends of the cylinder to form a torus with $K = L \times L$ demes.

In each generation, a given gene stays put with probability
$1-m$. Otherwise it migrates with equal probability to one of the
four demes to the immediate N, S, E, or W. The probability of
migrating to each of these demes is therefore $m/4$.

Why a torus? We are assuming that all migration is between adjacent
demes on a 2-dimensional surface. In other words, isolation by
distance in 2 dimensions. We could do this using a square lattice of
demes, but then things would get weird at the edges, because demes on
the edge of the lattice don't have as many neighbors. There is a
theorem (cite Seneta) showing that the effects of these edges become
less and less important as the number of demes grows. Consequently,
they probably don't matter in a real population of substantial
size. We can ignore these effects by wrapping the lattice around a
torus. 

    python3 msp.py -r 1>sim1.psmcfa 2>heterozygosity.txt
    python3 msp.py -r 1>sim2.psmcfa 2>>heterozygosity.txt

    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim1.psmc sim1.psmcfa
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim2.psmc sim2.psmcfa

    cat sim1.psmc sim2.psmc > torus.psmc

    ~/distrib/psmc/utils/psmc_plot.pl -g 29 -u 1.4e-8 -p \
	-M 'sim1,sim2' -f 'Helvetica,20' \
    -x 4e4 -X 4e6 -T "N=25, d=144, m=0.1" torus-144 torus.psmc

Reading the code in `psmc_plot.pl`, it appears that the mutation rate
specified by `-u` is per nucleotide site per generation. To verify
this, search Li's source code for `opts{u}`.



