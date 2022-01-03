#include "functions.h"
#include "MersenneTwister.h"
#include <cstdlib>
#include <iostream>
using namespace std;
using std::atoi;


// Random number generator:
MTRand rnd;



int main(int argc, char* argv[])
{
    

    // Parameters:
    
    int N, r, nP, nL, nbgen, nbgene, pas, sampleS, rep;
    double a, s, h, U, sP, hP, Un, L, linkP;
    
    
    N = atoi(argv[1]);
    a = atof(argv[2]);
    r = atoi(argv[3]);
    nL = atoi(argv[4]);
    sP = atof(argv[5]);
    hP = atof(argv[6]);
    linkP = atof(argv[7]);
    s = atof(argv[8]);
    h = atof(argv[9]);
    U = atof(argv[10]);
    L = atof(argv[11]);
    Un = atof(argv[12]);
    sampleS = atoi(argv[13]);
    nbgene = atoi(argv[14]);
    nbgen = atoi(argv[15]);
    pas = atoi(argv[16]);
    rep = atoi(argv[17]);
    
    // Simulation:
    
    recursion(N, a, r, nL, sP, hP, linkP, s, h, U, L, Un, sampleS, nbgene, nbgen, pas, rep);
    
    
}
