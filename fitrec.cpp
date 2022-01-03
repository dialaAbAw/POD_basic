#include "functions.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

extern MTRand rnd;

//fitness function calculates the multiplicative fitness

double fitness(vector<double> &c1, vector<double> &c2, double wHe, double wHo)
{
    int s1 = c1.size() - 1;
    int s2 = c2.size() - 1;
    
    double w = 1.0;
    
    while ((s1 > -1) || (s2 > -1))
    {

        if (s1 == -1) { w *= pow(wHe, s2 + 1); break;}
        else if (s2 == -1) { w *= pow(wHe, s1 + 1); break;}
        
        if(c1[s1] == c2[s2]){ w *= wHo; s1--; s2--;}
        else if(c1[s1] > c2[s2]){ w *= wHe; s1--;}
        else{ w *= wHe; s2--;}
        
    }

    return w;
    
}


// Measuring heterozygocity
int heter(vector<double> &c1, vector<double> &c2)
{
    int s1 = c1.size() - 1;
    int s2 = c2.size() - 1;
    int he = 0;
    
    while ((s1 > -1) || (s2 > -1))
    {
        
        if (s1 == -1) { he += s2 + 1; break;}
        else if (s2 == -1) { he += s1 + 1; break;}
        
        if(c1[s1] == c2[s2]){ s1--; s2--;}
        else if(c1[s1] > c2[s2]){ he++; s1--;}
        else{ he++; s2--;}
        
    }
    
    return he;
    
}

// rec function: generates recombinant genome segment "res"
// from parental genome segments "c1" and "c2"
// posneut contains the positions of neutral loci
// nbCo is the number of crossovers
void rec(chr &res, chr &c1, chr &c2, int nbCo)
{
	vector<double> pos;
	int j, locus;
    int nbsel1 = c1.sel.size();
    int nbsel2 = c2.sel.size();
    res.sel.clear();

    // The allele at the unlinked neutral locus is inherited from c1
    res.nlocus = c1.nlocus;
    
	// vector "pos" holds the positions of cross-overs:
  
    for (j = 0; j < nbCo; j++){ pos.push_back(rnd.rand());}
	sort(pos.begin(), pos.end());
    
    // loci under selection
    int locus1 = 0;
    int locus2 = 0;
    for (j = 0; j < nbCo; j++)
    {

        if (j % 2 == 0){
            while ((locus1 < nbsel1) && (c1.sel[locus1] < pos[j]))
            {
                res.sel.push_back(c1.sel[locus1]);
                locus1++;
            }
            while ((locus2 < nbsel2) && (c2.sel[locus2] < pos[j]))
            {
                locus2++;
            }
        }
        else{
            while ((locus2 < nbsel2) && (c2.sel[locus2] < pos[j]))
            {
                res.sel.push_back(c2.sel[locus2]);
                locus2++;
            }
            while ((locus1 < nbsel1) && (c1.sel[locus1] < pos[j]))
            {
                locus1++;
            }
            
        }
    }

    if (nbCo % 2 == 0)
        while (locus1 < nbsel1)
        {
            res.sel.push_back(c1.sel[locus1]);
            locus1++;
    
        }
    else
        while (locus2 < nbsel2)
        {
            res.sel.push_back(c2.sel[locus2]);
            locus2++;
        }

    sort(res.sel.begin(), res.sel.end());

    // Recombination in the POD zones
    nbsel1 = c1.pod.size();
    nbsel2 = c2.pod.size();
    res.pod.clear();

    locus1 = 0;
    locus2 = 0;
    for (j = 0; j < nbCo; j++)
    {

        if (j % 2 == 0)
        {
            while ((locus1 < nbsel1) && (c1.pod[locus1] < pos[j]))
            {
                res.pod.push_back(c1.pod[locus1]);
                locus1++;
            }
            while ((locus2 < nbsel2) && (c2.pod[locus2] < pos[j]))
            {
                locus2++;
            }
        }
        else
        {
            while ((locus2 < nbsel2) && (c2.pod[locus2] < pos[j]))
            {
                res.pod.push_back(c2.pod[locus2]);
                locus2++;
            }
            while ((locus1 < nbsel1) && (c1.pod[locus1] < pos[j]))
            {
                locus1++;
            }
        }
    }

    if (nbCo % 2 == 0)
        while (locus1 < nbsel1)
        {
            res.pod.push_back(c1.pod[locus1]);
            locus1++;
        }
    else
        while (locus2 < nbsel2)
        {
            res.pod.push_back(c2.pod[locus2]);
            locus2++;
        }

    sort(res.pod.begin(), res.pod.end());
}
