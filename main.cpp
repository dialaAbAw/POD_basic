#include "functions.h"
#include "MersenneTwister.h"
#include <cstdlib>
#include <iostream>
using namespace std;
using std::atoi;


// Random number generator:
MTRand rnd;

//int Nv, double av, int nLv, double sPv, double hPv, double linkPv, double Lv, int sampleSv, int nbgenv, int pasv, int repv, int Gv, double m


int main(int argc, char* argv[])
{


    // Parameters:

    int N, nL, nbgen, pas, sampleS, rep, G;
    double a, sP, hP, L, linkP, m;


    N = atoi(argv[1]);
    a = atof(argv[2]);
    //r = atoi(argv[3]);
    nL = atoi(argv[3]);
    //nEq = atof(argv[5]);
    sP = atof(argv[4]);
    hP = atof(argv[5]);
    linkP = atof(argv[6]);
    //s = atof(argv[9]);
    //h = atof(argv[10]);
    //U = atof(argv[11]);
    L = atof(argv[7]);
    //Un = atof(argv[13]);
    sampleS = atoi(argv[8]);
    //nbgene = atoi(argv[15]);
    nbgen = atoi(argv[9]);
    pas = atoi(argv[10]);
    rep = atoi(argv[11]);
    G = atoi(argv[12]);
    m = atof(argv[13]);

    // Simulation:

    recursion(N, a, nL, sP, hP, linkP, L, sampleS, nbgen, pas, rep, G, m);


}
