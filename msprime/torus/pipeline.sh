python3 msp.py -r 1>sim1.psmcfa 2>heterozygosity.txt
python3 msp.py -r 1>sim2.psmcfa 2>>heterozygosity.txt

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim1.psmc sim1.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o sim2.psmc sim2.psmcfa

cat sim1.psmc sim2.psmc > torus.psmc

~/distrib/psmc/utils/psmc_plot.pl -g 29 -u 1.4e-8 -p \
    -M 'sim1,sim2' -f 'Helvetica,20' \
    -x 2e4 -X 2e6 -Y 15 -T "2Nm = 2" torus torus.psmc
