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
 rv: 0 - Loci in the POD zone are equally spaced, 1 - Loci in the POD zone are randomly spaced
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
 */

void recursion(int Nv, double av, int rv, int nLv, double nEqv, double sPv, double hPv, double linkPv, double sv, double hv, double Uv, double Lv, double Unv, int sampleSv, int nbgenev, int nbgenv, int pasv, int repv)

{
    int i, j, k, nb, nb2, nb3, gen, mut, pa1, pa2, nbCo, nbFix, nbPFix, ns, repn;
    double rd, rd2, rd3, rd4, rdd1, rdd2, w, wbar, wmax, wfix, wp, wpod, div, nbdel, hebar, henb, wout, wself, wpout, wpself, hePbar, hePnb;

    int twoN = 2 * Nv;
    int twosampleS = 2*sampleSv;

    // The population is represented as a table containing haploid chromosomes
    chr pop[twoN];
    chr prepop[twoN];
    chr temp[twoN];
    chr ind1, ind2;

    // Vector to store individual identites to generate population susceptible to POD
    std::vector<int> Pos;

    // Heterozygote and Homozygote fitnesses of mutations involved in POD
    double Whetp = 1 - hPv * sPv;
    double Whomp = 1 - sPv;

    // Heterozygote and Homozygote fitnesses of deleterious mutations
    double Whet = 1 - hv * sv;
    double Whom = 1 - sv;

    // table of fitness values
    double Wij[Nv];

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
    /*std::vector<double> inda;
    std::vector<double> indb; */*
    std::vector<double> inds[Gv];   
    std::vector<double> indc;
    std::vector<double> fixedP;


for (repn = 1; repn <= repv; repn++){
 // Scaling physical distance in relation to the full chromosome map length
    rd = linkPv / (Lv);
// Placing the POD zone in the center of the chromosome, with nLv loci on either side of rd2 (a total of 2*nLv loci are positioned, nLv on each haplotype). This can be changed here.
    rd4 = 0.5;

    /*int nL = std::max(nLv, nLb);
    nLb = std::min(nLv, nLb); */
    

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
        {
            inds[g].clear();
        }

    if(rv == 0){
        //placing mutations at a distance 2*rd from one another on each haplotype
        for (j = 0; j < nL; j++){
            // Chromosomal region with POD zone overlap.
            fixedP.push_back(rd2);
            rd2 += rd;
            rd3 = (fixedP.size() - 1 ) %% Gv;     //rd3 = index de l'individu dans lequel ira la mutation. le -1 sert à retomber sur un rang de vecteur entre 0 et Gv-1 
            inds[rd3].push_back(rd2);
        }
    }
    
    else{
        // Random positioning of mutations on each haplotype 
        //total length of first POD zone
        rdd1 *=2;

        for (j = 0; j < nL; j++)
        {
            // placing mutations randomly between position rd2 and rd2 + rdd1

            do
            {
                rd4 = rnd.rand(rdd1);
                rd4 += rd2;
            } while (find(fixedP.begin(), fixedP.end(), rd4) != fixedP.end());

            fixedP.push_back(rd4);
            rd3 = (fixedP.size() - 1 ) %% Gv;
            inds[rd3].push_back(rd2);
         }
         
         for (g = 0; g < Gv; g++)
            {
                sort(inds[g].begin(), inds[g].end());
            }
    }

    fhout<< repn << " ";
    //First line of Haplo output file contains all the positions of each locus on each haplotype
       for(i = 0; i < nL; i++)
        fhout << inda[i] << " ";
    
       for(i = 0; i < nLb; i++)
        fhout << indb[i] << " ";
    
    fhout << endl;
    
    //Second line of Haplo output file indicates whether the mutation was initially on the first or second haplotype
    fhout << repn<<" ";

    for(i = 0; i < nL; i++)
        fhout << "1 ";

    for(i = 0; i < nLb; i++)
        fhout << "2 ";
    
    fhout << endl;
    
    
    

        // Re-initialising vectors with POD haplotype information
        fixedP.clear();


        //Burn-in simulations are run twice to obtain two different initial populations
        for (ns = 0; ns < 2; ns++)
        {

            ind1.sel.clear();
            ind1.pod.clear();
            ind1.nlocus = 0;
            fixedP.clear();
            nbFix = 0;
            for (j = 0; j < twoN; j++)
                pop[j] = ind1;

            // Recursion
            for (gen = 0; gen <= nbgenev; gen++)
            {
                // Introducing mutations
                for (i = 0; i < twoN; i++) // for each chromosome
                {
                    // number of new deleterious mutations
                    mut = poisdev(Uv);
                    for (j = 0; j < mut; j++)
                    {
                        // each mutation has a random position between 0 and 2L:
                        pop[i].sel.push_back(rnd.rand());
                    }
                    sort(pop[i].sel.begin(), pop[i].sel.end());

                    // At free recombining neutral locus (infinite alleles)
                    if (rnd.rand() < Unv)
                    {
                        pop[i].nlocus = rnd.rand();
                    }
                }

                //Measuring fitness
                wbar = 0;
                wmax = 0;
                for (j = 0; j < Nv; j++)
                {
                    nb = 2 * j;
                    // Fitness is calculated using "fitness" function defined in SelRec.cpp.
                    w = fitness(pop[nb].sel, pop[nb + 1].sel, Whet, Whom);

                    // Storing individual fitnesses
                    Wij[j] = w;

                    //Calculating mean population fitness
                    wbar += w;

                    // storing max fitness value
                    if (wmax < w)
                        wmax = w;
                }

                //Mean fitnessess taking fixed mutations into account
                wbar /= Nv;
                wbar *= pow(Whom, nbFix);

                //Normalising fitness values
                for (i = 0; i < Nv; i++)
                {
                    Wij[i] /= wmax;
                }

                //Reproduction
                for (j = 0; j < Nv; j++)
                {
                    nb = 2 * j;
                    // Selecting the maternal individual
                    do
                    {
                        i = rnd.randInt(Nv - 1);

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
                            i = rnd.randInt(Nv - 1);

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

                // Updating population
                for (i = 0; i < twoN; i++)
                    pop[i] = temp[i];

                //Removing fixed mutations every pasv generations
                if (gen % pasv == 0)
                {
                    for (i = pop[0].sel.size() - 1; i >= 0; i--)
                    {
                        j = 1;
                        rd = pop[0].sel[i];
                        while ((j < twoN) &&
                               (find(pop[j].sel.begin(), pop[j].sel.end(), rd) != pop[j].sel.end()))
                            j++;
                        if (j == twoN)
                        {
                            nbFix++;
                            fixedP.push_back(rd);
                            for (k = 0; k < twoN; k++)
                                pop[k].sel.erase(find(pop[k].sel.begin(), pop[k].sel.end(), rd));
                        }
                    }
                }
            }

            // Output of population statistics
            // Accounting for mutations fixed in pre-POD zones
            nb = nLv;
            wp = pow(Whomp, nb);
            wbar *= wp;
            wpod = wp;

            // Number of deleterious mutations per chromosome, compiling list of loci carrying mutations and calculating diversity at free recombing neutral locus
            Fr.clear();
            indc.clear();
            nbdel = 0;

            for (i = 0; i < twoN; i++)
            {

                k = pop[i].sel.size();
                nbdel += k;

                // Identifying all deleterious alleles segregating in the chromosome

                for (j = 0; j < k; j++)
                {
                    if (find(indc.begin(), indc.end(), pop[i].sel[j]) == indc.end())
                        indc.push_back(pop[i].sel[j]);
                }

                sort(indc.begin(), indc.end());

                //Identifying and counting neutral alleles
                nb = Fr.size();
                for (j = 0; j < nb; j++)
                    if (Fr[j].all == pop[i].nlocus)
                    {
                        Fr[j].freq++;
                        break;
                    }

                if (j == nb)
                {
                    alltemp.all = pop[i].nlocus;
                    alltemp.freq = 1;
                    Fr.push_back(alltemp);
                }
            }

            nb = Fr.size();
            div = 0;
            for (j = 0; j < nb; j++)
            {
                rd = Fr[j].freq / twoN;
                div += rd * rd;
            }

            nbdel /= double(twoN);

            // Number of mutations segregating in the population
            nb2 = indc.size();
            indc.clear();

            // Measuring heterozygocity in non POD zones
            hebar = 0;
            for (i = 0; i < Nv; i++)
            {
                j = 2 * i;
                rd = heter(pop[j].sel, pop[j + 1].sel);
                hebar += rd;
            }

            hebar /= Nv;
            henb = nb2;

            // Inbreeding depression calculated using sampleS individuals
            wout = 0;
            wself = 0;

            // Selfed offspring:
            for (j = 0; j < sampleSv; j++)
            {
                pa1 = 2 * (rnd.randInt(Nv - 1));

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

            wself /= sampleSv;

            // Outcrossed offspring:

            for (j = 0; j < sampleSv; j++)
            {
                pa1 = int(rnd.randInt(Nv - 1));

                do
                {
                    pa2 = rnd.randInt(Nv - 1);
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
                wout += w;
            }

            wout /= sampleSv;

            fout<< repn<<" " << ns - 2 << " " << wbar << " " << 1.0 - (wself / wout) << " " << nbFix << " " << nbdel << " " << hebar << " " << henb << " " << wpod << " " << 0 << " " << 0 << " " << 0 << " " << div << endl;

            // Ramdomly sampling individuals to create new population with POD zones included
            //Sampling from the first population (first simulation run, POD haplotype saved in vector inda)
            if (ns == 0)
            {
    
                nb = 0;
                // After admixture, each initial population would have contributed 50% of the individuals of the new population. This can be changed here (it must also be changed in the following "else" command)
                nb2 = 0.5 * double(Nv);
                indc = inda;
            }
            //Sampling from the second population (second simulation run, POD haplotype saved in vector indb)
            else
            {
                nb = 0.5 * double(Nv);
                nb2 = Nv;
                indc = indb;
            }

            pa1 = Nv - 1;
            // Re-initialising vector storing positions of sampled individuals that will be introduced into the admixed population
            Pos.clear();
            mut = fixedP.size();

            for (i = nb; i < nb2; i++)
            {
                j = 2 * i;
                do
                {
                    nb3 = rnd.randInt(pa1);
                } while (find(Pos.begin(), Pos.end(), nb3) != Pos.end());

                Pos.push_back(nb3);

                nb3 *= 2;

                prepop[j] = pop[nb3];
                for (k = 0; k < mut; k++)
                    prepop[j].sel.push_back(fixedP[k]);
                prepop[j].pod = indc;

                prepop[j + 1] = pop[nb3 + 1];
                for (k = 0; k < mut; k++)
                    prepop[j + 1].sel.push_back(fixedP[k]);
                prepop[j + 1].pod = indc;
            }

            indc.clear();
            fixedP.clear();
        }

        // Simulations with POD zones

        // Initialising population
        for (i = 0; i < twoN; i++)
            pop[i] = prepop[i];

        nbFix = 0;
        nbPFix = 0;
        fixedP.clear();

        // Recursion
        for (gen = 0; gen <= nbgenv; gen++)
        {
            //Removing fixed mutations and writing population statistics every pasv generations
            if (gen % pasv == 0)
            {

                // Fixed deleterious mutations
                for (i = pop[0].sel.size() - 1; i >= 0; i--)
                {
                    j = 1;
                    rd = pop[0].sel[i];
                    while ((j < twoN) &&
                           (find(pop[j].sel.begin(), pop[j].sel.end(), rd) != pop[j].sel.end()))
                        j++;
                    if (j == twoN)
                    {
                        nbFix++;
                        for (k = 0; k < twoN; k++)
                            pop[k].sel.erase(find(pop[k].sel.begin(), pop[k].sel.end(), rd));
                    }
                }

                // Fixed mutations in POD zones
                for (i = pop[0].pod.size() - 1; i >= 0; i--)
                {
                    j = 1;
                    rd = pop[0].pod[i];
                    while ((j < twoN) && (find(pop[j].pod.begin(), pop[j].pod.end(), rd) != pop[j].pod.end()))
                        j++;
                    if (j == twoN)
                    {
                       //Saving positions of mutations fixed in the POD zones
                        fixedP.push_back(rd);
                        for (k = 0; k < twoN; k++)
                            pop[k].pod.erase(find(pop[k].pod.begin(), pop[k].pod.end(), rd));
                    }
                }

                nbPFix = fixedP.size();
            }

            // Introducing mutations
            for (i = 0; i < twoN; i++) // for each chromosome
            {
                // number of new deleterious mutations
                mut = poisdev(Uv);
                for (j = 0; j < mut; j++)
                {
                    // each mutation has a random position between 0 and 2L:
                    pop[i].sel.push_back(rnd.rand());
                }
                sort(pop[i].sel.begin(), pop[i].sel.end());

                // At free recombining neutral locus (infinite alleles)
                if (rnd.rand() < Unv)
                {
                    pop[i].nlocus = rnd.rand();
                }
            }

            //Measuring fitness
            wbar = 0;
            wpod = 0;
            wmax = 0;
            for (j = 0; j < Nv; j++)
            {
                nb = 2 * j;
                // Fitness is calculated using "fitness" function defined in SelRec.cpp.
                w = fitness(pop[nb].sel, pop[nb + 1].sel, Whet, Whom);
                wp = fitness(pop[nb].pod, pop[nb + 1].pod, Whetp, Whomp);
                w *= wp;

                // Storing individual fitnesses
                Wij[j] = w;

                //Calculating mean population fitness
                wbar += w;
                wpod += wp;

                // storing max fitness value
                if (wmax < w)
                    wmax = w;
            }

            //Mean fitnessess taking fixed mutations into account
            wbar /= Nv;
            wbar *= pow(Whom, nbFix);
            wp = pow(Whomp, nbPFix);
            wbar *= wp;
            wpod /= Nv;
            wpod *= wp;

            //Normalising fitness values
            for (i = 0; i < Nv; i++)
                Wij[i] /= wmax;
            

            //Reproduction
            for (j = 0; j < Nv; j++)
            {
                nb = 2 * j;
                // Selecting the maternal individual
                do
                {
                    i = rnd.randInt(Nv - 1);

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
                        i = rnd.randInt(Nv - 1);

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


            //Removing fixed mutations and writing populaiton statistics every pasv generations
            if (gen % pasv == 0)
            {

                // Number of deleterious mutations per chromosome, compiling list of loci carrying mutations and calculating diversity at free recombing neutral locus
                Fr.clear();
                indc.clear();
                nbdel = 0;
                Pos.clear();

               //Randomly selecting sampleSv individuals on which to measure summary statistics
                for (i = 0; i < sampleSv; i++)
                {
                    do
                    {
                        nb3 = rnd.randInt(Nv - 1);
                        nb3 *= 2;

                    } while (find(Pos.begin(), Pos.end(), nb3) != Pos.end());

                    Pos.push_back(nb3);
                    Pos.push_back(nb3+1);
                }


                for (i = 0; i < twosampleS; i++)
                {
                    nb3 = Pos[i];
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
                }


                //Calculating neutral diversity

                nb = Fr.size();
                div = 0;
                for (j = 0; j < nb; j++)
                {
                    rd = Fr[j].freq / twosampleS;
                    div += rd * rd;
                }

                nbdel /= double(twosampleS);

                nb2 = indc.size(); 

                indc.clear();
               
                // number of polymorphic loci in POD zones
                for (i = 0; i < twosampleS; i++)
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
                indc.clear();

                // Measuring heterozygocity within the POD and throughout the genome
                hebar = 0;
                hePbar = 0;
                henb = 0;
                hePnb = 0;
                for (i = 0; i < sampleSv; i++)
                {
                    j = Pos[2*i];
                    rd = heter(pop[j].sel, pop[j + 1].sel);
                    rd2 = heter(pop[j].pod, pop[j + 1].pod);
                    hebar += rd;
                    hePbar += rd2;
                    
                }

                hebar /= sampleSv;
                henb = nb2 + nb; //toutes les mutations polymorphes dans la population
                hePbar /= sampleSv;
                hePnb = nb;

                          
               
                // Inbreeding depression calculated using sampleS individuals
                wout = 0;
                wself = 0;
                wpout = 0;
                wpself = 0;

                // Selfed offspring:
                for (j = 0; j < sampleSv; j++)
                {
                    pa1 = 2 * (rnd.randInt(Nv - 1));

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
                    wp = fitness(ind1.pod, ind2.pod, Whetp, Whomp);
                    w *= wp;
                    wself += w;
                    wpself += wp;
                }

                wself /= sampleSv;
                wpself /= sampleSv;

                // Outcrossed offspring:

                for (j = 0; j < sampleSv; j++)
                {
                    pa1 = int(rnd.randInt(Nv - 1));

                    do
                    {
                        pa2 = rnd.randInt(Nv - 1);
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
                    wp = fitness(ind1.pod, ind2.pod, Whetp, Whomp);
                    w *= wp;
                    wout += w;
                    wpout += wp;
                }

                wout /= sampleSv;
                wpout /= sampleSv;

                fout<< repn<<" " << gen << " " << wbar << " " << 1.0 - (wself / wout) << " " << nbFix << " " << nbdel << " " << hebar << " " << henb << " " << wpod << " " << 1.0 - (wpself / wpout) << " " << nbPFix << " " << hePbar << " " << hePnb << " " << div << endl;
                
            }

            // Updating population 
            for (i = 0; i < twoN; i++)
                pop[i] = temp[i];
        }

        //Writing state of mutations in POD zones in Haplo output file (fixed = 2, segregating = 1 and lost = 0)
        

        //Mutations initially present on first POD haplotype 
        fhout << repn<<" ";

        for (i = 0; i < nL; i++)
        {

            if ((fixedP.size() != 0) && (find(fixedP.begin(), fixedP.end(), inda[i]) != fixedP.end()))
            {
                fhout << "2 ";
                continue;
            }
            else
            {
                j = 0;
                while (j < twoN)
                {
                    if (find(pop[j].pod.begin(), pop[j].pod.end(), inda[i]) != pop[j].pod.end())
                    {
                        fhout << "1 ";
                        break;
                    }
                    j++;
                }

                if (j == twoN)
                {
                    fhout << "0 ";
                }
            }
        }

        //Mutations initially present on second POD haplotype 

        for (i = 0; i < nLb; i++)
        {

            if ((fixedP.size() != 0) && (find(fixedP.begin(), fixedP.end(), indb[i]) != fixedP.end()))
            {
                fhout << "2 ";
                continue;
            }
            else
            {
                j = 0;
                while (j < twoN)
                {
                    if (find(pop[j].pod.begin(), pop[j].pod.end(), indb[i]) != pop[j].pod.end())
                    {
                        fhout << "1 ";
                        break;
                    }  
                    j++;
                }

                if (j == twoN)
                {
                    fhout << "0 ";
                }
            }
        }

        fhout << endl;
}

}
