#ifndef PROCRYSTAL_NODE_H
#define PROCRYSTAL_NODE_H

#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"
#include "ring.h"
#include <random>

#include <iostream>
#include <SFML/Graphics.hpp>
#include <fstream>
#include <random>
#include <string>
#include <sstream>
#include <chrono>
#include <vector>
#include <sys/stat.h>
#include <iomanip>



using namespace std;

class Node {
    //Node in lattice with connectivity information and associated rings

private:

public:

    //Static variables
    static int autoId;

    //Data members
    int id;
    int environment,cluster;
    int numconnections;
    VecR<int> allowedCnxs,cnxs,rings;
    VecF<double> crd;
    int latticedimX, latticedimY;
    float fractionLlinear;
    //Constructors
    Node();
    Node(int maxCnxs, float fracl);

    //Member functions
    void resetRings();
    void setCrd(VecF<double> c);
    void addCnx(int cId);
    void addRing(int rId);
    void setOrientation(int orientation, int latdimX, int latdimY);
    void clockwiseCnxs(int latdimX, int latdimY);
    void randomCnxs(int numCnxs, mt19937 &mtGen, int latdimX, int latdimY);
    void randomCnxs(int numCnxs, VecF<int> pattern, mt19937 &mtGen);
    void bendorstraight(int dirs, int threecoordcount, int fourcoordcount, VecR<int> threecoord, VecR<int> fourcoord);

    void setCnxs(VecR<int> ids);
    void breakCnx(int id);
    void maximiseCnxs();
    void swapRing(int delId,int addId);
    VecR<int> getVacantCnxs();
    VecR<int> getReciprocalCnxs(VecR<Node> &others);
    VecR<int> getSharedRings(Node &other);
    VecR<int> getSharedCnxs(Ring &ring);


};


using namespace sf;


struct pipe{
public:
    static int autoId;
    VecR<Vector2i> dirs; //directions for a node
    VecR<int>cnxs;
    int id;
    Vector2i crds;
    int latticedimX, latticedimY;
    int orientation; //orientation of nodes
    float angle; bool on;


    pipe();
    pipe(int latdimX, int latdimY);
    void rotate();
    bool isConnect(Vector2i dir);

    void bendorstraighten();

    void setCnxs(VecR<int> Nodecnxs);
    void setOrientation();
    void setAngle();

    void cnxstodirs();

    string orientationstring();



//    bool checkneighbourdistance(Vector2i v, Vector2i d){
//        Vector2i neighbourpos = putin(v+d);
//        cout << neighbourpos.x - v.x <<","<<neighbourpos.y - v.y<< endl;
//        if (int(abs(neighbourpos.x - v.x))==1 or int(abs(neighbourpos.x - v.x))==N-1){
//            if (int(abs(neighbourpos.y - v.y))==1 or int(abs(neighbourpos.y - v.y))==N-1){
//                return true;
//            }
//            else return false;
//        }
//        else return false;
//    }

    int energynode(VecR<pipe> pipes, int latdimX, int latdimY);



};


//pipe grid[16][16];
//pipe& cell(Vector2i v) {return grid[v.x][v.y];}
//bool isOut(Vector2i v) {return !IntRect(0,0,16,16).contains(v);}
Vector2i putin(Vector2i neighbourpos, int latdimX, int latdimY);
int energytot(VecR<pipe> pipes, int latdimX, int latdimY);
int energyVector(Vector2i v, VecR<pipe> pipes, int latdimX, int latdimY);
bool isOut(Vector2i v, int latdimX, int latdimY);
vector<int> neighbours(int node, vector<int> display, int latdimX, int latdimY);
#endif //PROCRYSTAL_NODE_H