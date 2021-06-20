#ifndef PROCRYSTAL_SAMPLER_H
#define PROCRYSTAL_SAMPLER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <algorithm>

#include <sys/stat.h>

#include "outputfile.h"
#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"
#include "lattice.h"

class Sampler {
    //Generates procrystal samples

private:

    //Procrystal parameters
    string latticeType;
    bool pattern;
    VecF<int> patternR, patternL;
    int latticeDim, nodeCnd;
    float fractionH, fractionorderedH, fractionL, fractionorderedL, fractionLlinear;
    int usearray, useLlinear;
    string secondaryinputFile;
    int outStyle;
    string outputfolder;
    int mcswitch;
    int maxTry;
    //Monte Carlo parameters
    int randomSeed,numSamples;
    double temperature;

    //Output parameters
    string outputPrefix;
    bool writeSamples,writeRDFs,writeEnvs,writeSk,writeRingAn,writeChainAn;
    int skMaxN;
    double rdfDelta,skDelta;


public:

    //Constructor and setters
    Sampler();
    void setProcrystal(string latType, int cnd, bool pat, VecF<int> patR, VecF<int> patL, int latDim, float fracH, float fracordH, float fracL, float fracordL, int useL, float fracLlin, int usearr, string secinpt, string outfol);
    void setMonteCarlo(int seed, double temp, int samples, int mcsw, int maxTr);
    void setOutput(int outSty, string outPrefix, int write, int envs, int rdf, double delta,  int sk, double sdelta, int smaxn);

    //Start sampler
    void sample(Logfile &logfile);
    void bruteForce(Logfile &logfile);

};


#endif //PROCRYSTAL_SAMPLER_H
