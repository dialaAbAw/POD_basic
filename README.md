# POD_basic

C++ program for simulating pseudo-overdominance generated after an admixture event. Two initial dipoloid and partially self-fertilizing populations of equal size are simulated, each fixed for a different haplotype in the pseudo-ovordominant or POD zone. After a burn-in period, individuals from each population are sampled to generate a new admixed population. The simulation then continues to run for a defined number of generations, over which summary statistics are measured.

To compile this program the header file MersenneTwister.h is included. It is based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus, written by Richard J. Wagner.

The file main.cpp initialises the simulation and contains the list of input parameters, life_cycle.cpp contains the loop representing the life-cycle (processes of selection, measuring fitness, mutation and reproduction), the fitrec.cpp file contains the functions for measuring fitness, calculating heterozygosity and the recombination function, the header functions.h, contains definitions for different structures used and ranbin.cpp contains codes for distributions used within the program (Gaussian distribution, Poisson distribution, etc.).

To compile this program using GNU you need to type this command in the terminal, making sure you are in the right working directory (i.e. where all the necessary files are): g++ -o sims *.cpp *.h -lm

"sims" can be replaced by any other name you wish to name the executable.

Parameters are:
 N: Population size
 a: Selfing rate
 r: 0 - Loci in the POD zone are equally spaced, 1 - Loci in the POD zone are randomly spaced
 nL: Number of loci in each POD zone
 sP: selection coefficient of alleles in POD zone
 hP: dominance coefficient of alleles in POD zone
 linkP: Distance in Morgans between loci in POD zone
 U: deleterious mutation rate per haploid genome
 s: selection coefficient of mutations involved in Background Selection
 h: dominance coefficient of mutations involved in Background Selection
 L: genome map length (average number of cross-overs at meiosis)
 Un: Mutation rate for the unlinked neutral locus (infinite alleles)
 sampleS: sample size for estimating inbreeding depression and heterosis
 nbgene: Number of generations in pasv for burn-in time
 nbgen: Total number of generations
 pas: number of generations between output of simulation results
 rep: identifier number to avoid overwriting different simulations run with the same parameters
 

To launch, the command in the terminal is simply: ./sims N, a, r, nL, sP, hP, linkP, s, h, U, L, Un, sampleS, nbgene, nbgen, pas, rep
(replacing the letters by the value of the corresponding parameter) 
