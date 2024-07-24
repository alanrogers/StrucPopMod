# Circular stepping stone model 

Parameters:  $M = 4Nm = 10$, $N = 600$, $m = 5/4N$, $d=6$.

Steps in analysis:

    python msp.py -r 1>sim1.psmcfa 2>heterozygosity.txt
    python msp.py -r 1>sim2.psmcfa 2>>heterozygosity.txt
    python msp.py -r 1>sim3.psmcfa 2>>heterozygosity.txt

    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim1.psmc sim1.psmcfa &
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim2.psmc sim2.psmcfa &
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim3.psmc sim3.psmcfa

    psmcdata -u 1.4e-8 -g 29 sim1.psmc > sim1.txt
    psmcdata -u 1.4e-8 -g 29 sim2.psmc > sim2.txt
    psmcdata -u 1.4e-8 -g 29 sim3.psmc > sim3.txt

Then within julia:

    include("simplot.jl")
	savefig("simCircleM10.png")

