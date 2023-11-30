#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <vector>
#include <iostream>
#include "MersenneTwister.h"

using namespace std;



// definition of structure "chr" representing a chromosome:
// "sel" is a vector containing the positions of deleterious alleles along the chromosome
// "neut" is the table containing the identities at the nlocv neutral loci

struct chr
{
    vector<double> sel; // selected loci
    vector<double> pod; // loci involved in POD
    double nlocus; // neutral locus to calculate Ne and considered completely unlinked
};

// Storing allele identity and frequency
struct Nall
{
    double all;
    double freq;
};

// Prototypes of functions
void recursion(int Nv, double av, int nLv, double sPv, double hPv, double linkPv, double Lv, int sampleSv, int nbgenv, int pasv, int repv, int Gv, double m)
double gasdev();
double poisdev(const double xm);
double fitness(vector<double> &c1, vector<double> &c2, double wHe, double wHo);
void rec(chr &res, chr &c1, chr &c2, int nbco);
int heter(vector<double> &c1, vector<double> &c2);

#endif
