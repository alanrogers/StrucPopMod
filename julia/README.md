# imne
Effective population size under the island model of population structure

## Making graphs

The file txt/islandcircle.pdf was made as follows. First, cd into
src/StrucPopNe, then type "julia". Within the julia repl:

    ]activate .

The "]" gets you into Julia's "package mode". Then "activate ." makes
the StrucPopNe module accessible. After doing this, hit ^C to get back
into the Julia repl. Then type

    using Gadfly, Cairo, Fontconfig

Gadfly is below needed for "hstack". Cairo and Fontconfig are needed
for "PDF". Next:

For the islandcircle.pdf file:

    d = 6
    N = 600
	is = mkplot(:island, d, [6,7,8], 0, N,30;
	title="Island, d=$d N=$N X=0")
	ci = mkplot(:circle,d,[10,12,14],0,N,30;
	title="Circle, d=$d N=$N X=0")
    both = hstack(is, ci);
    draw(PDF("islandcircle.pdf", 30cm, 11cm), both)

Fst:

    julia> for M in [6,7,8]
               println(M, " ", fst(Val{:island}(), d, M, 0))
           end
    6 0.10373443983402489
    7 0.09025270758122744
    8 0.07987220447284345

    julia> for M in [10,12,14]
               println(M, " ", fst(Val{:circle}(), d, M, 0))
           end
    10 0.08860759493670886
    12 0.07494646680942184
    14 0.06493506493506493

Ratio of asymptotic Ne to initial Ne.

    julia> for M in [6,7,8]
               println(M, " ", ratio_of_change(Val{:island}(), d, M, 0))
           end
    6 6.709124377523067
    7 6.6061615472983295
    8 6.529277308114442

    julia> for M in [10,12,14]
               println(M, " ", ratio_of_change(Val{:circle}(), d, M, 0))
           end
    10 6.597580928441342
    12 6.49611583500809
    14 6.424075715529991

# Figures with extinction

Explored the case in which X=M, which implies that heterozysity w/i
groups is theta*(d+1)/2, roughly half the value in a model without
extinction. Setting this value = 2e-4, the observed value among
Neanderthals, gives N*(d+1) = 7142.857142857143. I have made several
graphs with values of N and d that satisfy this constraint.

Here is input for the island model:

    N=600
    600

    d=11
    11

    pisland = mkplot(:xisland, d, [4,5,6,7], 1,  N, 30;
    title="Island, d=$d, N=$N, X=M")

In this graph, M=X=6 generates an increase across 20ky, similar to
that in the Neanderthal data. For this value, 

    julia> fst(Val{:xisland}(), d, 6, 6)
    0.11210762331838565

    julia> ratio_of_change(Val{:xcircle}(), d, 6, 6)
    7.637422924807983

This model is plausible.
	
Same as above, but for circular stepping stone model.

    N = 600
	d = 11

    pcircle = mkplot(:xcircle, d, [10,20,30], 1,  N, 30;
    title="Circle, d=$d, N=$N, X=M")

    julia> fst(Val{:xcircle}(), d, 20, 20)
    0.06493506493506493
	
	ratio_of_change(Val{:xcircle}(), d, 20, 20)
    6.4693816470070145

Pretty close

    for M in [10,20,30]
	    println(M, " ", fst(Val{:xcircle}(), d, M, M))
	end
    10 0.12195121951219513
    20 0.06493506493506493
    30 0.04424778761061947

Plausible


Make a combined graph:

    both = hstack(pisland, pcircle);
    draw(PDF("xislandcircle.pdf", 30cm, 11cm), both)

Calculations for magnitude of change:

    julia> for M in [2,3,4,5,10,20,30]
               println(M, " ", ratio_of_change(xisland, 11, M, M))
           end
    2 8.653311931459035
    3 7.8921316273421445
    4 7.517200571083363
    5 7.294361719689457
    10 6.8541019662496865
    20 6.636997087674971
    30 6.56511990494422

    julia> for M in [2,3,4,5,10,20,30]
               println(M, " ", ratio_of_change(xcircle, 11, M, M))
           end
    2 11.312926484927337
    3 9.430410326521834
    4 8.519465958407334
    5 7.986233152643066
    10 6.958775250854595
    20 6.4693816470070145
    30 6.310544903560926

# Model in which X = 10M

    # mutation rate and heterozygosity
    u, h = 1.4e-8, 2e-4
    (1.4e-8, 0.0002)

    # X/M ratio
    a=10
    10

    # Find d for given values of other parameters
    get_d(h, u, 600, a)
    55.47619047619048

    N, d = 600, 55
    (600, 55)

    # make graph
    is10 = mkplot(:xisland, d, [4,6,8], a, N, 30;
	showNd=false,
    title="Island d=$d N=$N Nd=$(N*d) X=$(a)M")

    # M=6 looks best
    M = 6
    6

    fst(Val{:xisland}(), d, M, a*M)
	0.11967495690716573

    ratio_of_change(Val{:xisland}(), d, M, a*M)
    7.403044391478247

    N*d
    33000

    1/(a*M/(4N))
    40.0

    N, d, a = 600, 6, 0
	is0 = mkplot(:island, d, [4,6,8], 0, N,30;
	title="Island, d=$d N=$N Nd=$(N*d) X=0")

    # M=6 looks best
	M=6
	6

    fst(Val{:island}(), d, M, a*M)
    0.10373443983402489

    ratio_of_change(Val{:island}(), d, M, 0)
    6.709124377523067

    both = hstack(is0, is10);
    draw(PDF("island-lohiX.pdf", 30cm, 11cm), both)

