#ifndef PROCRYSTAL_LATTICE_H
#define PROCRYSTAL_LATTICE_H

#include <map>
#include <random>
#include <string>
#include <map>
#include <unordered_map>
#include <sstream>

#include <sys/stat.h>

#include <iostream>
#include <SFML/Graphics.hpp>
#include <fstream>
#include <vector>
#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"
#include "node.h"
#include "ring.h"
#include "chain.h"
#include "outputfile.h"

using namespace std;

class Lattice {

private:

    //Data members
    int latticeDimX, latticeDimY;
    int latticeCnd,nodeCnd;
    bool pattern;
    float fractionH, fractionorderedH, fractionL, fractionorderedL, fractionLlinear;
    int check4coord=0, check3coord=0, check2coord=0;
    VecF<int> patternR,patternL;
    VecR<Node> nodes;
    VecR<pipe> pipes;
    pipe& cell(Vector2i v){
        int ref = v.x+latticeDimX*v.y;
        return pipes[ref];
    }
    VecR<Ring> rings;
    VecR<Chain> chains;
    mt19937 mtGen;
    uniform_int_distribution<int> randRL;
    uniform_real_distribution<double> rand01;
    bool converged;
    VecF<double> periodicity;
    bool calcEnvs;
    string latticeCode;
    bool calcRings,calcChains;
    VecR<int> fourcoord;
    VecR<int> threecoord;
    VecR<int> twocoord;
    VecR<int> ncoords;
    VecR<int> ncoordsoriginal;
    int mcswitch;
    int useLlinear;

    int fourcoordcount;
    int threecoordcount;
    int twocoordcount;

    int lengthavaliablenodes;

    VecR<int> nodeenergy;

    int energy;
    int r2;

    int outStyle;
    int usearray;
    string secondaryinputFile;
    string outputfolder;
    string outputPrefix;

    //Member functions
    void initialiseSqLattice(int dimX, int dimY, bool restart=false);
    void initialiseTriLattice(int dim, bool restart=false);
    void initialiseSnubSqLattice(int dim, bool restart=false);
    void initialiseIsoSnubQuadLattice(int dim, bool restart=false);
    void initialiseTriHexLattice(int dim, bool restart=false);
    void initialiseHexLattice(int dim, bool restart=false);
    void threeandfourcoords();
    void initialiseProLattice();
    int findRings();
    int findChains();
    void findEnvironments();
    void calculateRDF(VecF<double> &x, VecF<double> &y, VecF<double> &rdf, double delta);
    void calculateSk(VecF<double> &x, VecF<double> &y, VecF<double> &sk, VecF<int> &multiplicity, double delta, int nMax);
    VecF<double> periodicCentreOfMass(VecF<double> &x, VecF<double> &y);

public:

    //Constructor
    Lattice();
    Lattice(int seed);

    //Member functions
    void initialise(string latType, int cnd, int latDimX, int latDimY, bool pat, float fracH, float fracordH, float fracL, float fracordL, int useL, float fracLlin, int usearr, string secinpt, string outfol, VecF<int> patR, VecF<int> patL, string prefix, int outsty, bool restart=false, bool crystal=false, bool envs=false);
    VecF<int> new_generate(int maxIterations, double temperature, int mcswitch);

    VecF<int> generate(int maxIterations, double temperature, int mcswitch);
    int generate(VecR<int> &pairA, VecR<int> &pairB);
    VecF<double> networkAnalysis(VecF<int> &k, VecF<int> &pk, VecF< VecF<int> > &ejk);
    void chainAnalysis(VecR<int> &chainLengths);
    VecF<double> calculateNetworkProperties(VecF<int> &k, VecF<int> &pk, VecF< VecF<int> > &ejk);
    void rdfAnalysis(VecF<double> &latticeRDF, VecF<double> &dualRDF, double delta);
    void nodeAnalysis(VecF<int> &pt, VecF< VecF<int> > &pst, VecF< VecF< VecF<int> > > &prst, OutputFile &testFile);
    void skAnalysis(VecF<double> &sk, VecF<int> &multiplicity, double delta, int nMax);
    int getNumNodes();
    int getNumActiveRings();
    VecF<double> getPeriodicity();
    VecF<int> getEnvironments();
    string getEnvironmentCode();
    void writeCrds(string prefix);
    void writeNetwork(string prefix, int index);
    void getEdgeCombinations(VecR<int> &pairA, VecR<int> &pairB, int &n, int &r);


    int optimumCnxs(int numCnxs);

    void generater2E();

    void checkCoords();

    int selectgoodnode();

    VecR<int> getavaliablenode();

    VecR<int> checkoptimumCnxs(int id);

    VecF<int> generateDisplay();


    int pipeenergytot(int latdim);
    int pipeenergyVector(Vector2i v, int latdim);
    int pipeenergynode(int pipeid, int latdim);

    void Monolayer(string prefix);
};


#endif //PROCRYSTAL_LATTICE_H
