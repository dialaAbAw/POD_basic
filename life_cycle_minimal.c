#include "functions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <cmath>
#include <csignal>
#include <cstring>
#include <algorithm>
#include <cstdio>
#include <iterator>
#include <unordered_set>
#include "MersenneTwister.h"

using namespace std;

extern MTRand rnd;

/* function cycle: iterates the life cycle until "equilibrium".
 Parameters are:
 Nv: Population size
 av: Selfing rate
 nLv: Number of deleterious mutations in the first POD zone
 nEqv: Proportion of nLv in second POD zone (i.e. nEqv = 1: both POD zones have the same number of deleterious mutations)
 sPv: selection coefficient of alleles in POD zone
 hPv: dominance coefficient of alleles in POD zone
 linkPv: Distance in Morgans beween loci in POD zone
 Uv: deleterious mutation rate per haploid genome
 sv: selection coefficient of mutations involved in Background Selection
 hv: dominance coefficient of mutations involved in Background Selection
 Lv: genome map length (average number of cross-overs at meiosis)
 Unv: Mutation rate for the unlinked neutral locus (infinite alleles)
 sampleSv: sample size for estimating inbreeding depression and heterosis
 nbgenev: Number of generations in pasv for burn-in time
 nbgenv: Total number of generations
 pasv: number of generations between output of simulation results
 repv: identifier number to avoid overwriting different simulations run with the same parameters
 Gv: demes index
 m: migration rate (probality for an offspring to descend from foreign parents)
 */

 
 //  A FAIRE : enlever var non appelées du call initial et des fout fhout
void recursion(int Nv, double av, int nLv, double nEqv, double sPv, double hPv, double linkPv, double sv, double hv, double Uv, double Lv, double Unv, int sampleSv, int nbgenev, int nbgenv, int pasv, int repv, int Gv, double m)

{
    int i, j, k, g, nb, nb2, nb3, gen, mut, pa1, pa2, nbCo, nbFix, nbPFix, ns, repn;
    double rd, rd2, rd3, rd4, rdd1, rdd2, w, wfix, wp, wpod, div, nbdel, hebar, henb, wout_intra, wout_inter, wself;


    int twoN = 2 * Nv;
    int twoNG = 2 * Nv *Gv;
    int NG =  Nv *Gv;
    int twosampleS = 2*sampleSv;

    // The population is represented as a table containing haploid chromosomes
    chr pop[twoNG];
    chr prepop[twoNG];
    chr temp[twoNG];
    chr ind1, ind2;

    // Vector to store individual identites to generate population susceptible to POD
    std::vector<int> Pos;

    // Heterozygote and Homozygote fitnesses of mutations involved in POD
    double Whetp = 1 - hPv * sPv;
    double Whomp = 1 - sPv;

    // table of fitness values
    double Wij[Nv];

    // table of mean fitness for each deme and for pod , and to store max fitness of each group
    double wbar[Gv];
    double wpod[Gv];
    double wmax[Gv];

    // Table to store neutral allele identities and frequencies at a given neutral locus
    std::vector<Nall> Fr;

    // Object used to store frequency of a neutral allele
    Nall alltemp;

    // creating Haplo output file storing positions and states of mutations in POD zones at the end of simulations (the name is created using parameter values, this can be modified here):
    char haplofile[256];
    stringstream nameH;
    nameH << "N" << Nv << "_a" << av << "_r" << rv <<"_nL" << nLv <<"_nEq"<< nEqv << "_link" << linkPv << "_sP" << sPv << "_hP" << hPv << "_U" << Uv << "_L" << Lv << "_Haplo.txt";
    nameH >> haplofile;
    ofstream fhout(haplofile);

    // creating output file containing summary statistics (the name is created using parameter values, this can be modified here):
    char mainfile[256];
    stringstream nameM;
     nameM << "N" << Nv << "_a" << av << "_r" << rv <<"_nL" << nLv <<"_nEq"<< nEqv << "_link" << linkPv << "_sP" << sPv << "_hP" << hPv << "_U" << Uv << "_L" << Lv << ".txt";
    nameM >> mainfile;
    ofstream fout(mainfile);

    // Creating POD zones: Positions of mutations are saved in vectors
    std::vector<double> inds[Gv]; //possibility to code this as a table
    std::vector<double> indc;
    std::vector<double> fixedP;



for (repn = 1; repn <= repv; repn++){
// Scaling physical distance in relation to the full chromosome map length
    rd = linkPv / (Lv);
// Placing the POD zone in the center of the chromosome, with nLv loci on either side of rd2 (a total of 2*nLv loci are positioned, nLv on each haplotype). This can be changed here.
    rd4 = 0.5;

    // Positioning mutations in the POD zone which is of map length 2*nLv*rd. For rv = 0 in the parameter set, mutations are equally spaced, if not, then they are placed randomly.
    // total number of initial deleterious mutations
    int nL = Gv*nLv
    // half length of first POD zone
    rdd1 = double(nL) *rd;
    // Starting position of POD on first haplotype
    rd2 = rd4 - rdd1;
    // End position of POD
    rd4 += rdd1;


    fixedP.clear();
    for (g=0; g<Gv; g++)
        inds[g].clear();


    //placing mutations at a distance 2*rd from one another on each haplotype
    for (j = 0; j < nL; j++){
        // Chromosomal region with POD zone overlap.
        fixedP.push_back(rd2);
        rd2 += rd;
        rd3 = (fixedP.size() - 1 ) % Gv;     //rd3 = index de l'individu dans lequel ira la mutation. le -1 sert à retomber sur un rang de vecteur entre 0 et Gv-1
        inds[rd3].push_back(rd2);
    }

    //First line of Haplo output file contains all the positions of each locus on each haplotype
    fhout<< repn << " ";
    for(g = 0; g < Gv; g++){
        for(i = 0; i < nLv; i++)
            fhout << inds[g][i] << " ";
    }
    fhout << endl;

    //Second line of Haplo output file indicates whether the mutation was initially on the first or second haplotype
    fhout << repn<<" ";

    for(g = 0; g < Gv; g++){
        for(i = 0; i < nLv; i++)
            fhout << g << " ";
    }
    fhout << endl;

    // Re-initialising vectors with POD haplotype information
    fixedP.clear();

    // assignation of each haplotype in each deme
    for (g=0; g<Gv; g++){
        j = g*twoN;                 
        for (i = j; i < j+twoN; i++){
            pop[i].pod = inds[g]
        }
    }


    // Recursion
    nbPFix = 0;
    for (gen = 0; gen <= nbgenv; gen++)
    {
        //Removing fixed mutations and writing population statistics every pasv generations
        if (gen % pasv == 0)
        {
            // Fixed mutations in POD zones
            for (i = pop[0].pod.size() - 1; i >= 0; i--)
            {
                j = 1;
                rd = pop[0].pod[i];
                while ((j < twoNG) && (find(pop[j].pod.begin(), pop[j].pod.end(), rd) != pop[j].pod.end()))
                    j++;
                if (j == twoNG)
                {
                   //Saving positions of mutations fixed in the POD zones
                    fixedP.push_back(rd);
                    for (k = 0; k < twoNG; k++)
                        pop[k].pod.erase(find(pop[k].pod.begin(), pop[k].pod.end(), rd));
                }
            }

            nbPFix = fixedP.size();
        }

        //Measuring fitness
        for (g = 0; g < Gv; g++)
        {
            for (j = 0; j < Nv; j++)
            {
                nb = 2 * j + g * twoN; // for each individuals of each deme
                // Fitness is calculated using "fitness" function defined in SelRec.cpp.
                w = fitness(pop[nb].pod, pop[nb + 1].pod, Whetp, Whomp);
                // Storing individual fitnesses
                nb = j + g * Nv
                Wij[nb] = w;  // nb is used as a placeholder here because it will be reinitialized in later loops
                wbar[g] += w;
                if (wmax[g] < w)
                    wmax[g] = w;
            }
            for (j = 0; j < Nv; j++)  // we need twice the same loop in order to use the real max value of wmax
            {
                nb = 2 * j + g * twoN;
                Wij[nb] /= wmax[g];
            }
            //Mean fitnessess taking fixed mutations into account
            wbar[g] /= Nv;
            wbar[g] *= pow(Whom, nbPFix);
        }

        //Reproduction
        for (g = 0; g < Gv; g++)
        {
            for (j = 0; j < Nv; j++)
            {
                nb = 2 * j + g * twoN;

                // ici je remplace g par k, et si il y a migration k est tiré différent de g
                if (rnd.rand() > m)
                {
                    k = g;  //no migration, the parents are selected in the same deme as the offspring
                }
                else
                {
                    do
                    {
                        k = rnd.randInt(Gv);
                    } while (k==g);
                }

                // Selecting the maternal individual
                do
                {
                    i = rnd.randInt(Nv - 1) + k * Nv ;

                } while (rnd.rand() > Wij[i]);
                pa1 = 2 * i;
                pa2 = 2 * i + 1;

                // sampling the number of crossovers
                nbCo = int(poisdev(Lv));

                // generating the first recombined chromosome
                if (0.5 < rnd.rand())
                    rec(ind1, pop[pa1], pop[pa2], nbCo);
                else
                    rec(ind1, pop[pa2], pop[pa1], nbCo);
                temp[nb] = ind1;


                // Selecting paternal individual

                // selfing:
                if (rnd.rand() < av)
                {
                    // sampling the number of crossovers
                    nbCo = int(poisdev(Lv));

                    if (0.5 < rnd.rand())
                        rec(ind1, pop[pa1], pop[pa2], nbCo);
                    else
                        rec(ind1, pop[pa2], pop[pa1], nbCo);
                    temp[nb + 1] = ind1;
                }

                // Outcrossing
                else
                {
                    do
                    {
                        i = rnd.randInt(Nv - 1) + Nv * k;

                    } while (rnd.rand() > Wij[i]);


                    pa1 = 2 * i;
                    pa2 = 2 * i + 1;

                    // sampling the number of crossovers
                    nbCo = int(poisdev(Lv));

                    // generating the first recombined chromosome
                    if (0.5 < rnd.rand())
                        rec(ind1, pop[pa1], pop[pa2], nbCo);
                    else
                        rec(ind1, pop[pa2], pop[pa1], nbCo);
                    temp[nb + 1] = ind1;
                }
            }
        }
    


    // writing population statistics every pasv generations
    if (gen % pasv == 0){

        // Number of deleterious mutations per chromosome, compiling list of loci carrying mutations and calculating diversity at free recombing neutral locus
        wout_intra = 0;
        wself = 0;
        wout_inter = 0; 

       //Randomly selecting sampleSv individuals in each deme on which to measure summary statistics
        for ( g = 0; g <Gv; g++)
        {
            // Inbreeding depression calculated using sampleS individuals

            // Selfed offspring:
            for (j = 0; j < sampleSv; j++)
            {
                pa1 = 2 * (rnd.randInt(Nv - 1))+g*Nv;

                //Recombination to create 1st chromosome
                // sampling the number of crossovers
                nbCo = int(poisdev(Lv));
                rd = rnd.rand();

                if (rd < 0.5)
                    rec(ind1, pop[pa1], pop[pa1 + 1], nbCo);
                else
                    rec(ind1, pop[pa1 + 1], pop[pa1], nbCo);

                //Recombination to create 2nd chromosome
                // sampling the number of crossovers
                nbCo = int(poisdev(Lv));

                rd = rnd.rand();
                if (rd < 0.5)
                    rec(ind2, pop[pa1], pop[pa1 + 1], nbCo);
                else
                    rec(ind2, pop[pa1 + 1], pop[pa1], nbCo);

                w = fitness(ind1.sel, ind2.sel, Whet, Whom);
                wself += w;
            }
            

            // Outcrossed offspring intra deme:

            for (j = 0; j < sampleSv; j++)
            {
                pa1 = int(rnd.randInt(Nv - 1))+g*Nv;

                do
                {
                    pa2 = rnd.randInt(Nv - 1)+g*Nv;
                } while (pa2 == pa1);

                pa1 *= 2;
                pa2 *= 2;

                // first chromosome
                // sampling the number of crossovers
                nbCo = int(poisdev(Lv));

                rd = rnd.rand();

                if (rd < 0.5)
                    rec(ind1, pop[pa1], pop[pa1 + 1], nbCo);
                else
                    rec(ind1, pop[pa1 + 1], pop[pa1], nbCo);

                // second chromosome
                // sampling the number of crossovers
                nbCo = int(poisdev(Lv));

                rd = rnd.rand();

                if (rd < 0.5)
                    rec(ind2, pop[pa2], pop[pa2 + 1], nbCo);
                else
                    rec(ind2, pop[pa2 + 1], pop[pa2], nbCo);

                w = fitness(ind1.sel, ind2.sel, Whet, Whom);
                wout_intra += w;
            }
            
            // Outcrossed offspring whole pop:

            for (j = 0; j < sampleSv; j++)
            {
                pa1 = int(rnd.randInt(Nv - 1))+g*Nv;

                do
                {
                    pa2 = rnd.randInt(NG - 1);
                } while (pa2 == pa1);

                pa1 *= 2;
                pa2 *= 2;

                // first chromosome
                // sampling the number of crossovers
                nbCo = int(poisdev(Lv));

                rd = rnd.rand();

                if (rd < 0.5)
                    rec(ind1, pop[pa1], pop[pa1 + 1], nbCo);
                else
                    rec(ind1, pop[pa1 + 1], pop[pa1], nbCo);

                // second chromosome
                // sampling the number of crossovers
                nbCo = int(poisdev(Lv));

                rd = rnd.rand();

                if (rd < 0.5)
                    rec(ind2, pop[pa2], pop[pa2 + 1], nbCo);
                else
                    rec(ind2, pop[pa2 + 1], pop[pa2], nbCo);

                w = fitness(ind1.sel, ind2.sel, Whet, Whom);
                wout_inter += w;
            }
        }
        
        wout_intra /= sampleSv * Gv; // repetition of sampleSv * Gv, maybe metre un placeholder name
        wout_inter /= sampleSv * Gv;
        wself /= sampleSv * Gv;
            
            Fr.clear();
            Pos.clear();
            indc.clear();
            nbdel = 0;
            for (i = 0; i < sampleSv* Gv; i++)
            {
                do
                {
                    nb3 = rnd.randInt(NG - 1);
                    nb3 *= 2;

                } while (find(Pos.begin(), Pos.end(), nb3) != Pos.end());

                Pos.push_back(nb3);
                Pos.push_back(nb3+1);

            }

            for (i = 0; i < twosampleS*Gv; i++) //meh 
            {
                nb = i ;
                nb3 = Pos[nb];
                k = pop[nb3].sel.size();
                nbdel += k;

              // number of polymorphic loci in the genome (excluding initial mutations introduced the POD zones)
                for (j = 0; j < k; j++)
                {
                    if (find(indc.begin(), indc.end(), pop[nb3].sel[j]) == indc.end())
                        indc.push_back(pop[nb3].sel[j]);
                }
                sort(indc.begin(), indc.end());  
                
                //Identifying and counting neutral alleles
                nb = Fr.size();
                for (j = 0; j < nb; j++)
                    if (Fr[j].all == pop[nb3].nlocus)
                    {
                        Fr[j].freq++;
                        break;
                    }
                if (j == nb)
                {
                    alltemp.all = pop[nb3].nlocus;
                    alltemp.freq = 1;
                    Fr.push_back(alltemp);              
                }
                
                nbdel /= double(twosampleS*Gv);
                nb2 = indc.size();
                indc.clear();

                // number of polymorphic loci in POD zones
                for (i = 0; i < twosampleS*Gv; i++)
                {
                    nb3 = Pos[i];
                    k = pop[nb3].pod.size();
                    for (j = 0; j < k; j++)
                    {
                        if (find(indc.begin(), indc.end(), pop[nb3].pod[j]) == indc.end())
                            indc.push_back(pop[nb3].pod[j]);
                    }

                    sort(indc.begin(), indc.end());
                }

                nb = indc.size();
                
                henb = nb2 + nb; 
            }
            fout << gen << " " << wself << " " << wout_intra << " " << wout_inter << " " << nbdel << " " << nbfix << " " << henb << endl;
        }
    // Updating population
            for (i = 0; i < twoN*Gv; i++)
                pop[i] = temp[i];
         
}
