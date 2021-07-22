#include "lattice.h"

sf::Vector2i SecondUp(0,-1);
sf::Vector2i SecondDown(0,1);
sf::Vector2i SecondRight(1,0);
sf::Vector2i SecondLeft(-1,0);
sf::Vector2i SecondDIR[4] = {SecondUp,SecondRight,SecondDown,SecondLeft};

Lattice::Lattice() {
    //Default constructor

    //Set random number generator
    mtGen.seed(0);
}


Lattice::Lattice(int seed) {
    //Constructor with random seed

    //Set random number generator
    mtGen.seed(seed);
    rand01=uniform_real_distribution<double>(0,1);
    randRL=uniform_int_distribution<int>(0,2);
}


void Lattice::initialise(string latType, int cnd, int latDimX, int latDimY, bool pat, float fracH, float fracordH, float fracL, float fracordL, int useL, float fracLlin, int usearr, string secinpt, string outfol, VecF<int> patR, VecF<int> patL, string prefix, int outSty,
                         bool restart, bool crystal, bool envs) {
    //Initialise lattice

    //Set variables
    ////////////////////
    fractionH=fracH;
    fractionorderedH = fracordH;
    fractionL=fracL;
    fractionorderedL = fracordL;

    useLlinear = useL;
    fractionLlinear = fracLlin;

    outStyle = outSty;
    usearray = usearr;
    secondaryinputFile = secinpt;
    outputPrefix = prefix;
    outputfolder = outfol;

    latticeDimX = latDimX;
    latticeDimY = latDimY;

    lengthavaliablenodes = 0;

    ////////////////////
    nodeCnd = cnd;

    pattern = pat;
    patternR = patR;
    patternL = patL;
    pipe::autoId = 0;
    Node::autoId = 0;
    Ring::autoId = 0;
    Ring::totalActive = 0;
    Chain::autoId = 0;
    converged=false;
    periodicity=VecF<double>(2);
    calcEnvs=envs;
    latticeCode = latType+to_string(nodeCnd);
    if(nodeCnd>2){
        calcRings=true;
        calcChains=false;
    }
    else{
        calcRings=false;
        calcChains=true;
    }
    //Make regular lattice
    if(latType=="sq") initialiseSqLattice(latDimX, latDimY,restart);
    else if(latType=="tri") initialiseTriLattice(latDimX, latDimY,restart);
    else if(latType=="snub") initialiseSnubSqLattice(latDimX, latDimY,restart);
    else if(latType=="isosnub") initialiseIsoSnubQuadLattice(latDimX, latDimY,restart);
    else if(latType=="trihex") initialiseTriHexLattice(latDimX, latDimY,restart);
    else if(latType=="hex") initialiseHexLattice(latDimX, latDimY,restart);

    //Calculate ring coordinates
    for(int i=0; i<rings.n; ++i){
        VecF<double> x(rings[i].nodes.n),y(rings[i].nodes.n);
        for(int j=0; j<rings[i].nodes.n; ++j){
            x[j] = nodes[rings[i].nodes[j]].crd[0];
            y[j] = nodes[rings[i].nodes[j]].crd[1];
        }
        VecF<double> com=periodicCentreOfMass(x,y);
        rings[i].setCrd(com);
    }

    //Write crystal
    if(!restart){
        for(int i=0; i<nodes.n; ++i) nodes[i].maximiseCnxs();
        writeCrds(prefix);
        writeNetwork(prefix,-1);
    }
    //Find Coordinations
    if(!crystal) threeandfourcoords();
    //Make initial procrystalline lattice
    if(!crystal) {
        cout << "INITIALISISE PRO LATTICE" << endl;
        initialiseProLattice();
    }
    //writeNetwork("./output/mc",-1);
}

void Lattice::initialiseSqLattice(int dim, bool restart) {
    //Initialise square lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 4;
    int dimSq = dim*dim;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq);
        pipes = VecR<pipe>(0,dimSq);

        for(int i=0; i<dimSq; ++i) {
            nodes.addValue(Node(latticeCnd, fractionLlinear));
            pipes.addValue(pipe(latticeDim));
            pipe &p = pipes[i];
            p.setCnxs(nodes[i].cnxs);
            p.cnxstodirs();
            p.setOrientation();
            p.setAngle();
        }
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,dimSq+(latticeCnd-nodeCnd)*dimSq);
    for(int i=0; i<dimSq; ++i) rings.addValue(Ring(4*dim));
    periodicity = dim;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id=0;
        for(int y=0; y<dim; ++y) {
            crd[1]=y;
            for(int x=0; x<dim; ++x) {
                crd[0]=x;
                nodes[id].setCrd(crd);
                ++id;
            }
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<dim; ++y){
            for(int x=0; x<dim; ++x){
                cnx = y * dim + (id + dim - 1) % dim;
                nodes[id].addCnx(cnx);
                cnx = (id + dim) % dimSq;
                nodes[id].addCnx(cnx);
                cnx = y * dim + (id + 1) % dim;
                nodes[id].addCnx(cnx);
                cnx = (id + dimSq - dim) % dimSq;
                nodes[id].addCnx(cnx);
                id += 1;
            }
        }
    }

    //Associate rings to nodes
    int id=0,ring;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dim; ++x){
            ring = id;
            nodes[id].addRing(ring);
            ring = y * dim + (id + 1) % dim;
            nodes[id].addRing(ring);
            ring=((y * dim + (id + 1) % dim) + dimSq - dim) % dimSq;
            nodes[id].addRing(ring);
            ring=(id + dimSq - dim) % dimSq;
            nodes[id].addRing(ring);
            id += 1;
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//BEWARE ! NODE() definition has been changed to allow FracL !! If you dont need this, remove it from node.cpp and remove the frac loop in randomCnxs
void Lattice::initialiseTriLattice(int dim, bool restart) {
    //Initialise triangular lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 6;
    int dimSq = dim*dim;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq);
        for(int i=0; i<dimSq; ++i) nodes.addValue(Node(latticeCnd, fractionLlinear));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,2*dimSq+(latticeCnd-nodeCnd)*2*dimSq);
    for(int i=0; i<2*dimSq; ++i) rings.addValue(Ring(2*4*dim));
    periodicity[0] = dim;
    periodicity[1] = dim*sqrt(3)/2;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id=0;
        double dy=sqrt(3.0)*0.5;
        for(int y=0; y<dim; ++y) {
            crd[1]=y*dy;
            for(int x=0; x<dim; ++x) {
                crd[0]=0.5*(y%2)+x;
                nodes[id].setCrd(crd);
                ++id;
            }
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<dim; ++y){
            for(int x=0; x<dim; ++x){
                cnx=y*dim+(id+dim-1)%dim;
                nodes[id].addCnx(cnx);
                if(y%2==0){
                    cnx=(cnx+dim)%dimSq;
                    nodes[id].addCnx(cnx);
                    cnx=(id+dim)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                else{
                    cnx=(id+dim)%dimSq;
                    nodes[id].addCnx(cnx);
                    cnx=((y+1)*dim+(id+1)%dim)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                cnx=y*dim+(id+1)%dim;
                nodes[id].addCnx(cnx);
                if(y%2==0){
                    cnx=(id+dimSq-dim)%dimSq;
                    nodes[id].addCnx(cnx);
                    cnx=(dimSq+(y-1)*dim+(id+dim-1)%dim)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                else{
                    cnx=(dimSq+(y-1)*dim+(id+dim+1)%dim)%dimSq;
                    nodes[id].addCnx(cnx);
                    cnx=(dimSq+(y-1)*dim+(id+dim)%dim)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                ++id;
            }
        }
    }

    //Associate rings to nodes
    int id=0,ring;
    int dimSq2=2*dimSq;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dim; ++x){
            ring=(2*y*dim+(id+dim-1)%dim);
            nodes[id].addRing(ring);
            ring=((2*y+1)*dim+(id)%dim);
            nodes[id].addRing(ring);
            ring=(2*y*dim+id%dim);
            nodes[id].addRing(ring);
            if(y%2==0){
                ring=(2*y*dim+id%dim-dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
                ring=(2*(y-1)*dim+(id+dim-1)%dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
                ring=(2*(y-1)*dim+(id+dim-1)%dim+dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
            }
            else{
                ring=(2*y*dim+(id+1)%dim-dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
                ring=(2*(y-1)*dim+id%dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
                ring=(2*(y-1)*dim+id%dim+dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
            }
            ++id;
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseSnubSqLattice(int dim, bool restart) {
    //Initialise snub-square lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 5;
    int xDim=dim,yDim=dim*2;
    int dimSq4=4*xDim*yDim;
    int dim4=4*xDim;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq4);
        for(int i=0; i<dimSq4; ++i) nodes.addValue(Node(latticeCnd, fractionLlinear));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,3*dimSq4/2+(latticeCnd-nodeCnd)*3*dimSq4/2);
    for(int i=0; i<3*dimSq4/2; ++i) rings.addValue(Ring(2*4*dim));
    double dxa=sqrt(3.0)/2;
    double dxe=1.0;
    double dya=0.5;
    double dye=0.5+sqrt(3.0)/2;
    periodicity[0] = xDim*(dxa*2+dxe);
    periodicity[1] = yDim*dye;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id=0;
        for(int y=0; y<yDim; ++y){
            crd[1]=dye*(y+0.5);
            for(int x=0; x<xDim; ++x){
                crd[0]=x*(2*dxa+dxe)+(y%2)*(dxe/2+dxa);
                nodes[y*dim4+x*4+0].setCrd(crd);
                crd[0]+=dxa;
                crd[1]+=dya;
                nodes[y*dim4+x*4+1].setCrd(crd);
                crd[1]-=2*dya;
                nodes[y*dim4+x*4+2].setCrd(crd);
                crd[0]+=dxa;
                crd[1]+=dya;
                nodes[y*dim4+x*4+3].setCrd(crd);
            }
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<yDim; ++y){
            for(int x=0; x<xDim; ++x){
                //environment 0
                cnx=y*dim4+(id+dim4-1)%dim4;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(y+1)*dim4+((id%dim4-2)+dim4)%dim4;
                else cnx=(id+dim4+2)%dimSq4;
                nodes[id].addCnx(cnx);
                cnx=id+1;
                nodes[id].addCnx(cnx);
                cnx=id+2;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(((y-1)*dim4+((id%dim4-3)+dim4)%dim4)+dimSq4)%dimSq4;
                else cnx=id-dim4+1;
                nodes[id].addCnx(cnx);
                ++id;
                //environment 1
                cnx=id-1;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(y+1)*dim4+((id%dim4-2)+dim4)%dim4;
                else cnx=(id+dim4+2)%dimSq4;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(y+1)*dim4+((id%dim4-1)+dim4)%dim4;
                else cnx=((y+1)*dim4+(id+3)%dim4)%dimSq4;
                nodes[id].addCnx(cnx);
                cnx=id+2;
                nodes[id].addCnx(cnx);
                cnx=id+1;
                nodes[id].addCnx(cnx);
                ++id;
                //environment 2
                cnx=id-2;
                nodes[id].addCnx(cnx);
                cnx=id-1;
                nodes[id].addCnx(cnx);
                cnx=id+1;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(((y-1)*dim4+((id%dim4-2)+dim4)%dim4)+dimSq4)%dimSq4;
                else cnx=(y-1)*dim4+((id%dim4+2)+dim4)%dim4;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(((y-1)*dim4+((id%dim4-3)+dim4)%dim4)+dimSq4)%dimSq4;
                else cnx=(y-1)*dim4+((id%dim4+1)+dim4)%dim4;
                nodes[id].addCnx(cnx);
                ++id;
                //environment 3
                cnx=id-1;
                nodes[id].addCnx(cnx);
                cnx=id-2;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(y+1)*dim4+((id%dim4-1)+dim4)%dim4;
                else cnx=((y+1)*dim4+((id%dim4+3)+dim4)%dim4)%dimSq4;
                nodes[id].addCnx(cnx);
                cnx=y*dim4+(id+dim4+1)%dim4;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(((y-1)*dim4+((id%dim4-2)+dim4)%dim4)+dimSq4)%dimSq4;
                else cnx=(y-1)*dim4+((id%dim4+2)+dim4)%dim4;
                nodes[id].addCnx(cnx);
                ++id;
            }
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate rings to nodes
    map<string,int> ringCodes;
    int id0,id1,id2,id3;
    int ringId=0;
    for(int i=0; i<nodes.n; ++i){
        id0=i;
        for(int j=0; j<5; ++j){
            id1=nodes[i].cnxs[j];
            id2=nodes[i].cnxs[(j+1)%5];
            VecR<int> ringPath(0,4);
            ringPath.addValue(i);
            ringPath.addValue(id1);
            ringPath.addValue(id2);
            if(!vContains(nodes[id1].cnxs,id2)){
                VecR<int> common=vCommonValues(nodes[id1].cnxs,nodes[id2].cnxs);
                common.delValue(id0);
                id3=common[0];
                ringPath.addValue(id3);
            }
            ringPath=vSort(ringPath);
            string rCode="";
            for(int j=0; j<ringPath.n; ++j) rCode += "#" + to_string(ringPath[j]);
            if(ringCodes.count(rCode)==0){
                ringCodes[rCode]=ringId;
                ++ringId;
            }
        }
    }
    for(int i=0; i<nodes.n; ++i) {
        id0 = i;
        for (int j = 0; j < 5; ++j) {
            id1 = nodes[i].cnxs[j];
            id2 = nodes[i].cnxs[(j + 1) % 5];
            VecR<int> ringPath(0, 4);
            ringPath.addValue(i);
            ringPath.addValue(id1);
            ringPath.addValue(id2);
            if (!vContains(nodes[id1].cnxs, id2)) {
                VecR<int> common = vCommonValues(nodes[id1].cnxs, nodes[id2].cnxs);
                common.delValue(id0);
                id3 = common[0];
                ringPath.addValue(id3);
            }
            ringPath = vSort(ringPath);
            string rCode= "";
            for(int j=0; j<ringPath.n; ++j) rCode += "#" + to_string(ringPath[j]);
            ringId=ringCodes.at(rCode);
            nodes[i].addRing(ringId);
        }
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseIsoSnubQuadLattice(int dim, bool restart) {
    //Initialise isosnub-quadrilateral lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 5;
    int dimSq=dim*dim;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq);
        for(int i=0; i<dimSq; ++i) nodes.addValue(Node(latticeCnd, fractionLlinear));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,3*dimSq/2+(latticeCnd-nodeCnd)*3*dimSq/2);
    for(int i=0; i<3*dimSq/2; ++i) rings.addValue(Ring(2*4*dim));
    periodicity[0] = dim;
    periodicity[1] = dim/2+dim*sqrt(3)/4;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id=0;
        for(int y=0; y<dim; ++y){
            for(int x=0; x<dim; ++x){
                crd[0]=x;
                if((y%4)>1) crd[0]+=0.5;
                nodes[id].setCrd(crd);
                ++id;
            }
            if(y%2==0) crd[1]+=1;
            else crd[1]+=sqrt(3)/2;
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<dim; ++y){
            for(int x=0; x<dim; ++x){
                cnx = y * dim + (id + dim - 1) % dim;
                nodes[id].addCnx(cnx);
                if(y%4==1){
                    cnx = ((id + dimSq + dim - 1) % dimSq);
                    if(x==0) cnx+=dim;
                    nodes[id].addCnx(cnx);
                }
                cnx = (id + dim) % dimSq;
                nodes[id].addCnx(cnx);
                if(y%4==3){
                    cnx = (id + dimSq + dim + 1) % dimSq;
                    if(x==dim-1) cnx-=dim;
                    nodes[id].addCnx(cnx);
                }
                cnx = y * dim + (id + 1) % dim;
                nodes[id].addCnx(cnx);
                if(y%4==2){
                    cnx = id - dim + 1;
                    if(x==dim-1) cnx-=dim;
                    nodes[id].addCnx(cnx);
                }
                cnx = (id + dimSq - dim) % dimSq;
                nodes[id].addCnx(cnx);
                if(y%4==0){
                    cnx = ((id + dimSq - dim - 1) % dimSq);
                    if(x==0) cnx+=dim;
                    nodes[id].addCnx(cnx);
                }
                ++id;
            }
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate rings to nodes
    map<string,int> ringCodes;
    int id0,id1,id2,id3;
    int ringId=0;
    for(int i=0; i<nodes.n; ++i){
        id0=i;
        for(int j=0; j<5; ++j){
            id1=nodes[i].cnxs[j];
            id2=nodes[i].cnxs[(j+1)%5];
            VecR<int> ringPath(0,4);
            ringPath.addValue(i);
            ringPath.addValue(id1);
            ringPath.addValue(id2);
            if(!vContains(nodes[id1].cnxs,id2)){
                VecR<int> common=vCommonValues(nodes[id1].cnxs,nodes[id2].cnxs);
                common.delValue(id0);
                id3=common[0];
                ringPath.addValue(id3);
            }
            ringPath=vSort(ringPath);
            string rCode="";
            for(int j=0; j<ringPath.n; ++j) rCode += "#" + to_string(ringPath[j]);
            if(ringCodes.count(rCode)==0){
                ringCodes[rCode]=ringId;
                ++ringId;
            }
        }
    }
    for(int i=0; i<nodes.n; ++i) {
        id0 = i;
        for (int j = 0; j < 5; ++j) {
            id1 = nodes[i].cnxs[j];
            id2 = nodes[i].cnxs[(j + 1) % 5];
            VecR<int> ringPath(0, 4);
            ringPath.addValue(i);
            ringPath.addValue(id1);
            ringPath.addValue(id2);
            if (!vContains(nodes[id1].cnxs, id2)) {
                VecR<int> common = vCommonValues(nodes[id1].cnxs, nodes[id2].cnxs);
                common.delValue(id0);
                id3 = common[0];
                ringPath.addValue(id3);
            }
            ringPath = vSort(ringPath);
            string rCode= "";
            for(int j=0; j<ringPath.n; ++j) rCode += "#" + to_string(ringPath[j]);
            ringId=ringCodes.at(rCode);
            nodes[i].addRing(ringId);
        }
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseTriHexLattice(int dim, bool restart) {
    //Initialise trihex lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 4;
    int dim3=3*dim;
    int dim_2=dim/2;
    int dimSq3_4=dim*dim3/4;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq3_4);
        for(int i=0; i<dimSq3_4; ++i) nodes.addValue(Node(latticeCnd, fractionLlinear));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,dimSq3_4+(latticeCnd-nodeCnd)*dimSq3_4);
    for(int i=0; i<dimSq3_4; ++i) rings.addValue(Ring(4*dim));
    periodicity[0] = dim;
    periodicity[1] = 2*sqrt(3)*dim/4;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        double dy=sqrt(3)/2;
        VecF<double> crd(2);
        int id=0;
        for(int y=0; y<dim/4; ++y){
            for(int yy=0; yy<4; ++yy) {
                if (yy % 4 == 0 || yy % 4 == 2) {
                    for (int x = 0; x < dim; ++x) {
                        crd[0] = x;
                        nodes[id].setCrd(crd);
                        ++id;
                    }
                }
                else if (yy % 4 == 1) {
                    for (int x = 0; x < dim_2; ++x) {
                        crd[0] = 2*x+0.5;
                        nodes[id].setCrd(crd);
                        ++id;
                    }
                }
                else if (yy % 4 == 3) {
                    for (int x = 0; x < dim_2; ++x) {
                        crd[0] = 2*x+1.5;
                        nodes[id].setCrd(crd);
                        ++id;
                    }
                }
                crd[1] += dy;
            }
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<dim/4; ++y){
            for(int yy=0; yy<4; ++yy){
                if(yy%4==0){
                    for(int x=0; x<dim; ++x){
                        cnx = y*dim3+(x+dim-1)%dim;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+int(floor(x/2.0));
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+(x+1)%dim;
                        nodes[id].addCnx(cnx);
                        cnx = (y*dim3-dim_2+(dim+int(ceil(x/2.0))-1)%dim_2+dimSq3_4)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        ++id;
                    }
                }
                if(yy%4==1){
                    for(int x=0; x<dim_2; ++x){
                        cnx = y*dim3+dim+dim_2+2*x;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+dim_2+2*x+1;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+2*x+1;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+2*x;
                        nodes[id].addCnx(cnx);
                        ++id;
                    }
                }
                if(yy%4==2){
                    for(int x=0; x<dim; ++x){
                        cnx = y*dim3+dim+dim_2+(x+dim-1)%dim;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+dim_2+dim+(int(ceil(x/2.0))-1+dim_2)%dim_2;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+dim_2+(x+1)%dim;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+int(floor(x/2.0));
                        nodes[id].addCnx(cnx);
                        ++id;
                    }
                }
                if(yy%4==3){
                    for(int x=0; x<dim_2; ++x){
                        cnx = ((y+1)*dim3+2*x+1)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        cnx = ((y+1)*dim3+(2*x+2+dim)%dim)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        cnx = ((y+1)*dim3-dim_2-dim+(2*x+2+dim)%dim+dimSq3_4)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        cnx = ((y+1)*dim3-dim_2-dim+(2*x+1+dim)%dim+dimSq3_4)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        ++id;
                    }
                }
            }
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate rings to nodes
    int id=0;
    int ring;
    for(int y=0; y<dim/4; ++y){
        for(int yy=0; yy<4; ++yy){
            if(yy%4==0){
                for(int x=0; x<dim; ++x){
                    ring = y*dim3+(x+dim-1)%dim;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+(int(ceil(x/2.0)))%dim_2;
                    nodes[id].addRing(ring);
                    ring = y*dim3+(x+0)%dim;
                    nodes[id].addRing(ring);
                    ring = (y*dim3-dim_2+int(floor(x/2.0))+dimSq3_4)%dimSq3_4;
                    nodes[id].addRing(ring);
                    ++id;
                }
            }
            if(yy%4==1){
                for(int x=0; x<dim_2; ++x){
                    ring = y*dim3+dim+x;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+dim_2+2*x;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+(dim_2+x+1)%dim_2;
                    nodes[id].addRing(ring);
                    ring = y*dim3+2*x;
                    nodes[id].addRing(ring);
                    ++id;
                }
            }
            if(yy%4==2){
                for(int x=0; x<dim; ++x){
                    ring = y*dim3+dim+dim_2+(x+dim-1)%dim;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+dim_2+dim+int(floor(x/2.0));
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+dim_2+(x+0)%dim;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+(int(ceil(x/2.0)))%dim_2;
                    nodes[id].addRing(ring);
                    ++id;
                }
            }
            if(yy%4==3){
                for(int x=0; x<dim_2; ++x){
                    ring = (y+1)*dim3-dim_2+x;
                    nodes[id].addRing(ring);
                    ring = ((y+1)*dim3+(2*x+1+dim)%dim)%dimSq3_4;
                    nodes[id].addRing(ring);
                    ring = (y+1)*dim3-dim_2+(dim_2+x+1)%dim_2;
                    nodes[id].addRing(ring);
                    ring = ((y+1)*dim3-dim_2-dim+(2*x+1+dim)%dim+dimSq3_4)%dimSq3_4;
                    nodes[id].addRing(ring);
                    ++id;
                }
            }
        }
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseHexLattice(int dim, bool restart) {
    //Initialise hexagonal lattice

    //Initialise nodes and rings
    latticeCnd = 3;
    int dimX = dim/2;
    int dimSq = dim*dimX;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq);
        for(int i=0; i<dimSq; ++i) nodes.addValue(Node(latticeCnd, fractionLlinear));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,2*dimSq+(latticeCnd-nodeCnd)*2*dimSq);
    for(int i=0; i<dimSq/2; ++i) rings.addValue(Ring(dim));
    periodicity[0] = dimX*sqrt(3);
    periodicity[1] = 3*dim/4;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id = 0;
        double dx = sqrt(3.0);
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dimX; ++x) {
                crd[0] += dx;
                nodes[id].setCrd(crd);
                ++id;
            }
            if (y % 4 == 0) {
                crd[0] = dx / 2;
                crd[1] += 0.5;
            } else if (y % 4 == 1) {
                crd[0] = dx / 2;
                crd[1] += 1.0;
            } else if (y % 4 == 2) {
                crd[0] = 0;
                crd[1] += 0.5;
            } else {
                crd[0] = 0;
                crd[1] += 1.0;
            }
        }

        //Make node connections
        id = 0;
        int cnx;
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dimX; ++x) {
                if (y % 4 == 0) {
                    cnx = (y + 1) * dimX + (x + dimX - 1) % dimX;
                    nodes[id].addCnx(cnx);
                    cnx = (y+1)*dimX+x;
                    nodes[id].addCnx(cnx);
                    cnx=(dimSq+(y-1)*dimX+x)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                else if (y % 4 == 1) {
                    cnx = (y - 1) * dimX + x;
                    nodes[id].addCnx(cnx);
                    cnx = (y - 1) * dimX + (x + 1) % dimX;
                    nodes[id].addCnx(cnx);
                    cnx = (y + 1) * dimX + x;
                    nodes[id].addCnx(cnx);
                }
                else if (y % 4 == 2) {
                    cnx = (y + 1) * dimX + x;
                    nodes[id].addCnx(cnx);
                    cnx = (y + 1) * dimX + (x + 1) % dimX;
                    nodes[id].addCnx(cnx);
                    cnx = (y - 1) * dimX + x;
                    nodes[id].addCnx(cnx);
                }
                else {
                    cnx = (y - 1) * dimX + (x + dimX - 1) % dimX;
                    nodes[id].addCnx(cnx);
                    cnx = (y-1)*dimX+x;
                    nodes[id].addCnx(cnx);
                    cnx=((y+1)*dimX+x)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                ++id;
            }
        }
    }

    //Associate rings to nodes
    int id=0,ring;
    int dimR=dimX;
    int dimRSq=dimR*dimR;
    int set4=0;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dimX; ++x){
            if(y%4==0){
                ring=set4*dimR*2+x;
                nodes[id].addRing(ring);
                ring=((set4+dimR/2-1)%(dimR/2))*dimR*2+dimR+(x+dimR-1)%dimR;
                nodes[id].addRing(ring);
                ring=((set4+dimR/2-1)%(dimR/2))*dimR*2+dimR+x;
                nodes[id].addRing(ring);
            }
            else if(y%4==1){
                ring=set4*dimR*2+x;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+(x+1)%dimR;
                nodes[id].addRing(ring);
                ring=((set4+dimR/2-1)%(dimR/2))*dimR*2+dimR+x;
                nodes[id].addRing(ring);
            }
            else if(y%4==2){
                ring=set4*dimR*2+dimR+x;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+x;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+(x+1)%dimR;
                nodes[id].addRing(ring);
            }
            else{
                ring=set4*dimR*2+dimR+(x+dimR-1)%dimR;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+dimR+x;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+x;
                nodes[id].addRing(ring);
            }
            ++id;
        }
        if(y%4==3) ++set4;
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}

///////////////////////////////////////////////////////////////////////

void Lattice::threeandfourcoords() {

    if (usearray == 1) {
        //cout << secondaryinputFile << endl;
        ifstream secondinputFile("./"+secondaryinputFile+".dat", ios::in);
        if(!secondinputFile.good()) cout << "Cannot find input file procrystal.inpt" << endl;
        string skip,line;
        int numberlines=0;
        fourcoordcount = int(nodes.n * fractionH);
        twocoordcount = int(nodes.n * fractionL);
        threecoordcount = nodes.n - fourcoordcount - twocoordcount;
        int newfourcoordcount = 0;
        int newthreecoordcount = 0;
        int newtwocoordcount = 0;
        int dummy;
        while (getline(secondinputFile, line)){
            ++numberlines;
            istringstream(line)>>dummy;
            //cout << dummy << endl;
            dummy = int(dummy);
            if (dummy == 2){
                newtwocoordcount += 1;
            }
            else if (dummy == 4){
                newfourcoordcount += 1;
            }
            else if (dummy == 30 or dummy == 31 or dummy ==3){
                newthreecoordcount += 1;
            }
            else{
                cout << "UNKNOWN VALUE" << endl;
            }
        }
        cout << "               N+1 :  " << newfourcoordcount  << "   :::    " << fourcoordcount << endl;
        cout << "               N   :  " << newthreecoordcount << "   :::    " << threecoordcount << endl;
        cout << "               N-2 :  " << newtwocoordcount   << "   :::    " << twocoordcount << endl;
        cout << "               --------------------------" << endl;
        cout << "               TOT :  " << fourcoordcount + threecoordcount + twocoordcount << endl;
        secondinputFile.close();



//        cout << "CHECK FILE LENGTH" << endl;
        if (latticeDim*latticeDim != numberlines){
            cout << "INPUT FILE HAS DIFFERENT DIMENSIONS TO EXPECTED" << endl;
            latticeDim = sqrt(numberlines);
            cout << "DIMENSIONS UPDATED but THIS WILL BREAK due to node properties" << endl;
        }
        else{
//            cout << "We're good to go" << endl;
        }
        cout << endl;
        cout << "CHECKING THE RATIOS" << endl;
        if (newfourcoordcount != fourcoordcount or newthreecoordcount != threecoordcount or newtwocoordcount != twocoordcount){
            cout << "INPUT FILE HAS DIFFERENT 4/3/2 RATIO TO EXPECTED" << endl;
            fourcoordcount = newfourcoordcount;
            threecoordcount = newthreecoordcount;
            twocoordcount = newtwocoordcount;
            cout << "RATIO  UPDATED - this should be ok ..." << endl;
        }
        else{
//            cout << "RATIOS ARE ACTUALLY OK" << endl;
        }

        ifstream thirdinputFile("./"+secondaryinputFile+".dat", ios::in);

        fourcoord = VecR<int>(fourcoordcount);
        threecoord = VecR<int>(threecoordcount);
        twocoord = VecR<int>(twocoordcount);
        ncoords = VecR<int>(nodes.n);
        ncoordsoriginal = VecR<int>(nodes.n);

        int k2=0;
        int k3=0;
        int k4=0;
        int i=0;
        while (getline(thirdinputFile, line)){

            istringstream(line)>>dummy;
            dummy = int(dummy);
            //cout << dummy << endl;
            if (dummy == 2){
                twocoord[k2] = i;
                k2+=1;
            }
            else if (dummy == 4){
                fourcoord[k4] = i;
                k4+=1;
            }
            else if (dummy == 30 or dummy == 31 or dummy ==3){
                threecoord[k3] = i;
                k3+=1;
            }
            else{
                cout << "UNKNOWN VALUE" << endl;
            }

            if (dummy == 30 or dummy == 31 or dummy == 3){
                ncoords[i] = 3;
                ncoordsoriginal[i] = 3;
            }
            else{
                ncoords[i] = dummy;
                ncoordsoriginal[i] = dummy;
            }
            i++;
        }

        cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;

        cout << "               FRACTION N+1 :  " << fractionH << endl;
        cout << "               FRACTION N   :  " << 1 - fractionH - fractionL << endl;
        cout << "               FRACTION N-1 :  " << fractionL << endl;
        cout << endl;
        cout << "               FRACTION ORDERED N+1  : " << fractionorderedH << endl;
        cout << "               FRACTION ORDERED N-1  : " << fractionorderedL << endl;
        cout << endl;
        cout << "               FRACTION LINEAR    :  " << fractionLlinear << endl;
        cout << endl;



    }

    else {

        cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;

        cout << "               FRACTION N+1 :  " << fractionH << endl;
        cout << "               FRACTION N   :  " << 1 - fractionH - fractionL << endl;
        cout << "               FRACTION N-1 :  " << fractionL << endl;
        cout << endl;
        cout << "               FRACTION ORDERED N+1  : " << fractionorderedH << endl;
        cout << "               FRACTION ORDERED N-1  : " << fractionorderedL << endl;
        cout << endl;
        cout << "               FRACTION LINEAR    :  " << fractionLlinear << endl;
        cout << endl;

        fourcoordcount = int(nodes.n * fractionH);
        twocoordcount = int(nodes.n * fractionL);
//    threecoordcount = int(nodes.n * (1-fractionH-fractionL));
        threecoordcount = nodes.n - fourcoordcount - twocoordcount;

        cout << fourcoordcount << "    " << threecoordcount << "    " << twocoordcount << "    "
             << nodes.n - fourcoordcount - threecoordcount - twocoordcount << endl;

        if (twocoordcount % 2 != 0) {
            twocoordcount += 1;
            threecoordcount -= 1;
        }
        if (fourcoordcount % 2 != 0) {
            fourcoordcount += 1;
            threecoordcount -= 1;
        }
        if ((nodes.n - fourcoordcount - threecoordcount - twocoordcount) == -1) {
            threecoordcount -= 1;
        } else if ((nodes.n - fourcoordcount - threecoordcount - twocoordcount) == +1) {
            threecoordcount += 1;
        }

        cout << fourcoordcount << "    " << threecoordcount << "    " << twocoordcount << "    "
             << nodes.n - fourcoordcount - threecoordcount - twocoordcount << endl;


        int dummy;
        int runningcount = 0;
        int orderedfourcoord, orderedthreecoord, orderedtwocoord;
        bool compare;


        fourcoord = VecR<int>(fourcoordcount);
        threecoord = VecR<int>(threecoordcount);
        twocoord = VecR<int>(twocoordcount);
        ncoords = VecR<int>(nodes.n);
        ncoordsoriginal = VecR<int>(nodes.n);


        orderedfourcoord = fourcoordcount * fractionorderedH;
        orderedthreecoord = threecoordcount * fractionorderedL;

        cout << "               N+1 :  " << fourcoordcount << endl;
        cout << "               N   :  " << threecoordcount << endl;
        cout << "               N-2 :  " << twocoordcount << endl;
        cout << "               --------------------------" << endl;
        cout << "               TOT :  " << fourcoordcount + threecoordcount + twocoordcount << endl;
        int check4coord = 0, check3coord = 0, check2coord = 0;
        for (int i = 0; i < nodes.n; i++) {
            int coordinateno = nodes[i].cnxs.n;
            if (coordinateno == 4) {
                check4coord += 1;
            } else if (coordinateno == 3) {
                check3coord += 1;
            } else {
                check2coord += 1;
            }
        }
        cout << "CHECK 2,3,4 COORDINATE    " << check4coord << "    " << check3coord << "    " << check2coord << "    "
             << nodes.n - check4coord - check3coord - check2coord << endl;

        cout << endl;

////    cout << nodeCnd << endl;
////    cout << int(nodeCnd + 1) << endl;
////    cout << int(nodeCnd - 1) << endl;
////    cout << endl;
        //uniform_int_distribution<int> randNode(orderedfourcoord,nodes.n-1);
        uniform_int_distribution<int> randNode(orderedfourcoord, nodes.n - 1 - orderedthreecoord);
        //////////////////////////////
        //// if you want the 3 coord sites ordered at the other end, swap for the next two lines
        //orderedthreecoord= threecoordcount*fractionordered;
        //uniform_int_distribution<int> randNode(orderedfourcoord,nodes.n-1-orderedthreecoord);
        //////////////////////////////

        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
        ///  adding three coord list ///
        /////////////////////////////////////////////////////////////////////////////////////////
        int i = 0;
        while (i < fourcoordcount) {

            if (i < orderedfourcoord) {
                fourcoord[i] = i;
                ncoords[i] = int(nodeCnd + 1);
                ncoordsoriginal[i] = int(nodeCnd + 1);
                i += 1;
            } else {
                compare = true;
                dummy = randNode(mtGen);
//            cout << dummy << endl;
                for (int j = 0; j < i; j++) {
                    if (dummy == fourcoord[j]) {
                        compare = false;
                    }
                }
                if (compare) {
                    fourcoord[i] = dummy;
                    ncoords[dummy] = int(nodeCnd + 1);
                    ncoordsoriginal[dummy] = int(nodeCnd + 1);
                    i += 1;
                }
            }
        }

////////////////    cout << "After 4 coords" << endl;
////////////////    cout << ncoords << endl;

        int k;
//    i=nodes.n-1;
        i = 0;
        while (i < threecoordcount) {
            //while (i > nodes.n - threecoordcount - twocoordcount - 1) {

            k = nodes.n - 1 - i;

////////////////        cout << i << "/" << nodes.n - threecoordcount - 1 << endl;
////////////////        cout << i << "               " << nodes.n - orderedthreecoord -1 << "         " << nodes.n - i - 1<< endl;

            if (k > nodes.n - orderedthreecoord - 1) {
//            cout << "ORDERED 3" << endl;
                threecoord[i] = k;
                ncoords[k] = int(nodeCnd);
                ncoordsoriginal[k] = int(nodeCnd);
                i += 1;
            } else {
////            cout << "through ordered" << endl;
                compare = true;
                dummy = randNode(mtGen);
                // this is slow step ...
                for (int j = 0; j < fourcoordcount; j++) {
                    // if proposed atom is already 4 coord
                    if (dummy == fourcoord[j]) {
                        compare = false;
                    }
                }
////            cout << "past four coord sites" << endl;
                for (int j = 0; j < i; j++) {
                    // if proposed atom is already in the 3 coord list
                    if (dummy == threecoord[j]) {
                        compare = false;
                    }
                }
////            cout << compare << endl;
                if (compare) {
                    threecoord[i] = dummy;
                    ncoords[dummy] = int(nodeCnd);
                    ncoordsoriginal[dummy] = int(nodeCnd);
                    i += 1;
                }
            }
        }

////////////////////////    cout << "After 3 coords" << endl;
////////////////////////    cout << ncoords << endl;

        // we've sorted out the ordered species now
        //two coordinate sites sit between the two


        i = 0;
        while (i < twocoordcount) {

            compare = true;
            dummy = randNode(mtGen);
//        cout << dummy << endl;
            for (int j = 0; j < fourcoordcount; j++) {
                // search through disordered 4 coordinate sites
                if (dummy == fourcoord[j]) {
                    compare = false;
                }
            }
            for (int j = 0; j < threecoordcount; j++) {
                // search through disordered 3 coordinate sites
                if (dummy == threecoord[j]) {
                    compare = false;
                }
            }
            for (int j = 0; j < i; j++) {
                if (dummy == twocoord[j]) {
                    compare = false;
                }
            }

            if (compare) {
                twocoord[i] = dummy;
                ncoords[dummy] = int(nodeCnd - 1);
                ncoordsoriginal[dummy] = int(nodeCnd - 1);
                i += 1;
            }
        }
    }


////    cout << "After 2 coords" << endl;
////    cout << ncoords << endl;

//    for (i=0;i<nodes.n-1;i++){
//        cout << i << "    " << ncoords[i] << endl;
//    }

}



/////////////////// NON RANDOM EQUIVALENT ////////////////////////////////////

//    VecR<int> fourcoord (nodes.n);
//    VecR<int> threecoord (nodes.n);
//
//
//    int fourcoordcount = 0;
//    int threecoordcount = 0;
//
//    uniform_int_distribution<int> randNode(0,nodes.n-1);
//
//    for (int i = 0; i < nodes.n; ++i) {
////        if (i%2==0) {
////        if (i < 0.5 * nodes.n) {
////            if (i!=20 and i!=21 and i!=40 and i!=41){
////            if (20>=i or i>44){
//        if (i == 20 or i == 21){
//            fourcoord[fourcoordcount] = i;
//            fourcoordcount +=1;
//        }
//        else {
//            threecoord[threecoordcount] = i;
//            threecoordcount += 1;
//        }
//    }


/////////////////////////////////////////////////////////////////////////////

void Lattice::checkCoords(){
    check4coord=0, check3coord=0, check2coord=0;
    for (int i=0;i<nodes.n;i++) {
        int coordinateno = nodes[i].cnxs.n;
        if (coordinateno == 4) {
            check4coord += 1;
        }
        else if (coordinateno == 3) {
            check3coord += 1;
        }
        else if (coordinateno == 2) {
            check2coord += 1;
        }
    }
}

void Lattice::initialiseProLattice() {
    //Remove random connections from crystal lattice to make procrystalline

    cout << "INITIALISING PRO LATTICE" << endl;

    //Write starting crystalline structure
    for(int i=0; i<nodes.n; ++i) nodes[i].maximiseCnxs();

    cout << endl;

    checkCoords();
    cout <<  "CHECK 2,3,4 COORDINATE    " << check4coord << "    " << check3coord << "    " << check2coord << "    " << nodes.n - check4coord - check3coord - check2coord << endl;

    for (int i=0;i<twocoordcount;i++){
        if (vContains(threecoord, twocoord[i])){
            cout << "Three coord contains two coord value" << endl;
        }
        else if (vContains(twocoord, threecoord[i])){
            cout << "Two coord contains three coord value" << endl;
        }
    }

    cout << endl;

    if(!pattern) {
        for (int i = 0; i < fourcoordcount; i++) {
            nodes[fourcoord[i]].randomCnxs(nodeCnd + 1, mtGen, latticeDim);
        }
        for (int i = 0; i < threecoordcount; i++) {
            nodes[threecoord[i]].randomCnxs(nodeCnd, mtGen, latticeDim);
        }
        for (int i = 0; i < twocoordcount; i++) {
            nodes[twocoord[i]].randomCnxs(nodeCnd - 1, mtGen, latticeDim);
        }
    }
    else{
        if(randRL(mtGen)==0){
            for (int i=0; i<fourcoordcount; i++) {
                nodes[fourcoord[i]].randomCnxs(nodeCnd + 1, patternR, mtGen);
            }
            for (int i=0; i<threecoordcount; i++) {
                nodes[threecoord[i]].randomCnxs(nodeCnd, patternR, mtGen);
            }
            for (int i=0; i<twocoordcount; i++) {
                nodes[twocoord[i]].randomCnxs(nodeCnd - 1, patternR, mtGen);
            }
        }
        else{
            for (int i=0; i<fourcoordcount; i++) {
                nodes[fourcoord[i]].randomCnxs(nodeCnd + 1, patternL, mtGen);
            }
            for (int i=0; i<threecoordcount; i++) {
                nodes[threecoord[i]].randomCnxs(nodeCnd, patternL, mtGen);
            }
            for (int i=0; i<twocoordcount; i++) {
                nodes[twocoord[i]].randomCnxs(nodeCnd - 1, patternL, mtGen);
            }
        }
    }
    cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    cout << "UPDATING PIPES" << endl;

    //
    //reflect these changes in the pipes
    for (int i=0;i<nodes.n;++i){
        int gridx,gridy;
        gridx = int(i%latticeDim);
        gridy = int((i-gridx)/latticeDim);
        pipe &p = pipes[i];
        p.id = nodes[i].id;

        if ((gridx+gridy)%2==0){
            nodes[i].setOrientation(1, latticeDim);
        }
        else {
            nodes[i].setOrientation(3,latticeDim);
        }

        p.setCnxs(nodes[i].cnxs);

//        cout << "Node Cnxs count = " << nodes[i].cnxs.n << endl;
//        for (int z=0;z<nodes[i].cnxs.n;z++){
//            cout << nodes[i].cnxs[z] << "  ";
//        }
//        cout << endl;
//
//        cout << "Node Cnxs count = " << pipes[i].cnxs.n << endl;
//        for (int z=0;z<pipes[i].cnxs.n;z++){
//            cout << pipes[i].cnxs[z] << "  ";
//        }
//        cout << endl;

        p.cnxstodirs();
        p.setOrientation();
        p.setAngle();
        p.crds = Vector2i {gridx,gridy};

    }

    cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    //

    for (int i=0;i<threecoordcount;i++){
        for (int j=0;j<twocoordcount;j++){
            if (threecoord[i] == twocoord[j]){
                cout << "Duplicate " << i << endl;
            }
        }
    }

    checkCoords();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "line 1540 " << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    cout << "               N+1 :  " << fourcoordcount << "    :    " << check4coord <<endl;
    cout << "               N   :  " << threecoordcount << "    :    " << check3coord << endl;
    cout << "               N-2 :  " << twocoordcount << "    :    " << check2coord << endl;
    cout << "               --------------------------" << endl;
    cout << "               TOT :  " << fourcoordcount + threecoordcount + twocoordcount << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}

void Lattice::generater2E() {
//    cout << "GENERATE E AND R2" << endl;
    nodeenergy = VecR<int>(nodes.n);
    energy=0;
    r2=0;
    VecR<int> rvec (2);
    for (int i=0; i<nodes.n; i++){
        nodeenergy[i] = ncoords[i]-nodes[i].getReciprocalCnxs(nodes).n;
    }
    for (int i=0; i<nodes.n; i++){
        energy +=ncoords[i]-nodes[i].getReciprocalCnxs(nodes).n;
        for (int j=(i+1); j<nodes.n; j++){
            if(nodeenergy[i]!=0 && nodeenergy[j]!=0){
                for (int k=0; k<2; k++){
                    rvec[k] = nodes[i].crd[k] - nodes[j].crd[k];
                    if(rvec[k]>latticeDim*0.5) rvec[k]-=latticeDim;
                    else if(rvec[k]<(-latticeDim*0.5)) rvec[k] +=latticeDim;
                    r2 += rvec[k]*rvec[k];
                }
            }
        }
    }
//    cout << "DONE" << endl;
}

int Lattice::optimumCnxs(int id) {
    // energetically/r2 best configuration for id and surrounding 4 nodes
    int nmax = 4;
    int numCnxs = ncoords[id];

//    cout << "OPTIMUM CONNECTIONS" << endl;
//    cout << "ID : " << id << endl;
    int delE=0;
    if (numCnxs == 4) return 0;
    if (numCnxs == 3) {
        VecR<int> energyVec(4);
        VecR<int> r2Vec(4);

//        cout << "NUM CNXS = 3" << endl;

        for (int i = 0; i < 4; i++) {
            // 4 orientations
            VecR<int> dummycnxs(4);
            dummycnxs = nodes[id].allowedCnxs;
            dummycnxs.delValue(nodes[id].allowedCnxs[i]);
            nodes[id].setCnxs(dummycnxs);

            generater2E();
            energyVec[i] = energy;
            r2Vec[i] = r2;
        }
//        cout << "ENERGY VEC" << endl;
//        cout << energyVec << endl;
//        cout << endl;
//        cout << "Min val" << endl;
//        cout << r2Vec.minimumValue() << endl;

        int to_remove = r2Vec.minimumValue()[0];
//        cout << "DELETE : " << to_remove << endl;
        VecR<int> dummycnxs (4);
        dummycnxs = nodes[id].allowedCnxs;
        dummycnxs.delValue(nodes[id].allowedCnxs[to_remove]);
        nodes[id].setCnxs(dummycnxs);
        delE = energyVec[to_remove];
    }
    if (numCnxs == 2) {
//        cout << "NUM CNXS = 2" << endl;
        VecR<int> energyVec(6);
        VecR<int> r2Vec(6);
        for (int i=0; i<4; i++) {
            for (int j = i + 1; j < 4; j++) {
                int ref = (4-i)*i + j-i-1;
                VecR<int> dummycnxs(4);
                dummycnxs = nodes[id].allowedCnxs;
                dummycnxs.delValue(nodes[id].allowedCnxs[i]);
                dummycnxs.delValue(nodes[id].allowedCnxs[j]);
                nodes[id].setCnxs(dummycnxs);

                generater2E();
                energyVec[ref] = energy;
                r2Vec[ref] = r2;
            }
        }
//        cout << "ENERGY VEC" << endl;
//        cout << energyVec << endl;
//        cout << endl;
//        cout << "Min val" << endl;
//        cout << r2Vec.minimumValue() << endl;
        int to_remove = r2Vec.minimumValue()[0];
//        cout << "DELETE : " << to_remove << endl;

        VecR<int> dummycnxs (4);
        dummycnxs = nodes[id].allowedCnxs;
        if(to_remove==0){
            dummycnxs.delValue(nodes[id].allowedCnxs[0]);
            dummycnxs.delValue(nodes[id].allowedCnxs[1]);
        }
        else if(to_remove==1){
            dummycnxs.delValue(nodes[id].allowedCnxs[0]);
            dummycnxs.delValue(nodes[id].allowedCnxs[2]);
        }
        else if(to_remove==2){
            dummycnxs.delValue(nodes[id].allowedCnxs[0]);
            dummycnxs.delValue(nodes[id].allowedCnxs[3]);
        }
        else if(to_remove==3){
            dummycnxs.delValue(nodes[id].allowedCnxs[1]);
            dummycnxs.delValue(nodes[id].allowedCnxs[2]);
        }
        else if(to_remove==4){
            dummycnxs.delValue(nodes[id].allowedCnxs[1]);
            dummycnxs.delValue(nodes[id].allowedCnxs[3]);
        }
        else if(to_remove==5){
            dummycnxs.delValue(nodes[id].allowedCnxs[2]);
            dummycnxs.delValue(nodes[id].allowedCnxs[3]);
        }
        nodes[id].setCnxs(dummycnxs);
        delE = energyVec[to_remove];
    }
    return delE;
}

VecR<int> Lattice::checkoptimumCnxs(int id) {
    // energetically/r2 best configuration for id and surrounding 4 nodes
    int nmax = 4;
    int numCnxs = ncoords[id];
    VecR<int> returnvals (2);
//    cout << "OPTIMUM CONNECTIONS" << endl;
//    cout << "ID : " << id << endl;
    int delE=0;
    if (numCnxs == 4) return 0;
    if (numCnxs == 3) {
        VecR<int> energyVec(4);
        VecR<int> r2Vec(4);

//        cout << "NUM CNXS = 3" << endl;

        for (int i = 0; i < 4; i++) {
            // 4 orientations
            VecR<int> dummycnxs(4);
            dummycnxs = nodes[id].allowedCnxs;
            dummycnxs.delValue(nodes[id].allowedCnxs[i]);
            nodes[id].setCnxs(dummycnxs);

            generater2E();
            energyVec[i] = energy;
            r2Vec[i] = r2;
        }
        int to_remove = r2Vec.minimumValue()[0];
        returnvals[0] = energyVec[to_remove];
        returnvals[1] = r2Vec[to_remove];
    }
    if (numCnxs == 2) {
        VecR<int> energyVec(6);
        VecR<int> r2Vec(6);
        for (int i=0; i<4; i++) {
            for (int j = i + 1; j < 4; j++) {
                int ref = (4-i)*i + j-i-1;
                VecR<int> dummycnxs(4);
                dummycnxs = nodes[id].allowedCnxs;
                dummycnxs.delValue(nodes[id].allowedCnxs[i]);
                dummycnxs.delValue(nodes[id].allowedCnxs[j]);
                nodes[id].setCnxs(dummycnxs);

                generater2E();
                energyVec[ref] = energy;
                r2Vec[ref] = r2;
            }
        }

        int to_remove = r2Vec.minimumValue()[0];
        returnvals[0] = energyVec[to_remove];
        returnvals[1] = r2Vec[to_remove];
    }
    return returnvals;
}

int Lattice::selectgoodnode(){
    VecR<int> avaliablenodes (nodes.n);
    // determine useful bond switches
    lengthavaliablenodes=nodes.n;

    for (int i=0;i<nodes.n;i++){
        avaliablenodes[i] = i;
    }

    //find nodes neighbouring unsaturated connections
    for (int i=0; i<nodes.n; i++){
        if(nodeenergy[i]==0){
            VecR<int> neighbours = nodes[i].allowedCnxs;
            bool check= true;
            for(int j=0;j<4;j++){
                if(nodeenergy[neighbours[j]] != 0){
                    check=false;
                }
            }
            if (check){
                avaliablenodes.delValue(i);
                lengthavaliablenodes -=1;
            }
        }
    }

    uniform_int_distribution<int> AvaliablerandNode(0,lengthavaliablenodes-1);

    int id = avaliablenodes[AvaliablerandNode(mtGen)];
    return id;
}

VecR<int> Lattice::getavaliablenode(){
    VecR<int> avaliablenodes (nodes.n);
    // determine useful bond switches
    lengthavaliablenodes=nodes.n;

    for (int i=0;i<nodes.n;i++){
        avaliablenodes[i] = i;
    }

    //find nodes neighbouring unsaturated connections
    for (int i=0; i<nodes.n; i++){
        if(nodeenergy[i]==0){
            VecR<int> neighbours = nodes[i].allowedCnxs;
            bool check= true;
            for(int j=0;j<4;j++){
                if(nodeenergy[neighbours[j]] != 0){
                    check=false;
                }
            }
            if (check){
                avaliablenodes.delValue(i);
                lengthavaliablenodes -=1;
            }
        }
    }
    return avaliablenodes;
}




VecF<int> Lattice::generateDisplay(){

    cout << "BEGINNING GENERATE DISPLAY" << endl;
    int N = latticeDim;
    double scale = double(double(16)/double(N));
    cout << "SCALE : " << scale << endl;
    std::vector<Vector2i> nodesCrds;
    int nCnxs;
//    int numberlines = 256;
    int numberlines = N*N;

    for (int i=0;i<nodes.n;i++) {
        nCnxs = nodes[i].numconnections;
        int x = int((numberlines) % N);
        int y = (int((numberlines - x) / N));
//        nodesCrds.push_back(Vector2i(x, y));
        pipe &p = pipes[i];

        p.setCnxs(nodes[i].cnxs);
        p.cnxstodirs();
        p.setOrientation();
        p.setAngle();
//        nodes[i].cnxs = nodes[i].allowedCnxs;
//        dummy1 = nodes[i].cnxs[0] - i;
//        dummy2 = nodes[i].cnxs[1] - i;
//        dummy3 = nodes[i].cnxs[2] - i;
//        dummy4 = nodes[i].cnxs[3] - i;
//        dummysum = (dummy1+dummy2+dummy3+dummy4)/latticedim;
//
//        VecR<int>dummy(4);
//        dummy[0]=dummy1;
//        dummy[1]=dummy2;
//        dummy[2]=dummy3;
//        dummy[3]=dummy4;
//
//        if (dummysum==0){
//            for (int j=0;j<4;j++){
//               if (dummy[j]==-1){
//                   A=j;
//               }
//               else if (dummy[j]==latticedim){
//                   D=j;
//               }
//               else if (dummy[j]==1){
//                   B=j;
//               }
//               else if (dummy[j]==-latticedim){
//                   C=j;
//               }
//               else {
//                   cout << "FAILS" << endl;
//               }
//            }
//        }
//        else if (dummysum > 0) {
//            // 1,2,3,6
//            if (dummysum == 1) {
//    //            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 6 "  << endl;
//                //A
//                for (int j=0;j<4;j++){
//                    if (dummy[j]==latticedim-1){
//                        A=j;
//                    }
//                    else if (dummy[j]==-latticedim){
//                        D=j;
//                    }
//                    else if (dummy[j]==1){
//                        B=j;
//                    }
//                    else if (dummy[j]==latticedim){
//                        C=j;
//                    }
//                    else {
//                        cout << "FAILS" << endl;
//                    }
//                }
//            }
//            else if (dummysum == latticedim-1){
//    //            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 1 "  << endl;
//                for (int j=0;j<4;j++){
//                    if (dummy[j]==-1){
//                        A=j;
//                    }
//                    else if (dummy[j]==latticedim){
//                        D=j;
//                    }
//                    else if (dummy[j]==1-latticedim){
//                        B=j;
//                    }
//                    else if (dummy[j]==latticedim*(latticedim-1)){
//                        C=j;
//                    }
//                    else {
//                        cout << "FAILS" << endl;
//                        cout << dummy[j] << endl;
//                    }
//                }
//            }
//            else if (dummysum == latticedim){
//    //            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 2 "  << endl;
//                for (int j=0;j<4;j++){
//                    if (dummy[j]==-1){
//                        A=j;
//                    }
//                    else if (dummy[j]==latticedim){
//                        D=j;
//                    }
//                    else if (dummy[j]==1){
//                        B=j;
//                    }
//                    else if (dummy[j]==latticedim*(latticedim-1)){
//                        C=j;
//                    }
//                    else {
//                        cout << "FAILS" << endl;
//                    }
//                }
//
//            }
//            else if (dummysum == latticedim+1){
//    //            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 3 "  << endl;
//                for (int j=0;j<4;j++){
//                    if (dummy[j]==latticedim-1){
//                        A=j;
//                    }
//                    else if (dummy[j]==latticedim){
//                        D=j;
//                    }
//                    else if (dummy[j]==1){
//                        B=j;
//                    }
//                    else if (dummy[j]==latticedim*(latticedim-1)){
//                        C=j;
//                    }
//                    else {
//                        cout << "FAILS" << endl;
//                    }
//                }
//            }
//        }
//        else {
//            // 5,7,8,9
//            if (dummysum == -1) {
//    //            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 5 "  << endl;
//                for (int j=0;j<4;j++){
//                    if (dummy[j]==-1){
//                        A=j;
//                    }
//                    else if (dummy[j]==latticedim){
//                        D=j;
//                    }
//                    else if (dummy[j]==1-latticedim){
//                        B=j;
//                    }
//                    else if (dummy[j]==-latticedim){
//                        C=j;
//                    }
//                    else {
//                        cout << "FAILS" << endl;
//                    }
//                }
//            }
//            else if (dummysum == 1-latticedim){
//    //            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 9 "  << endl;
//                for (int j=0;j<4;j++){
//                    if (dummy[j]==latticedim-1){
//                        A=j;
//                    }
//                    else if (dummy[j]==-latticedim*(latticedim-1)){
//                        D=j;
//                    }
//                    else if (dummy[j]==1){
//                        B=j;
//                    }
//                    else if (dummy[j]==-latticedim){
//                        C=j;
//                    }
//                    else {
//                        cout << "FAILS" << endl;
//                        cout << dummy[j] << endl;
//                    }
//                }
//            }
//            else if (dummysum == -latticedim){
//    //            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 8 "  << endl;
//                for (int j=0;j<4;j++){
//                    if (dummy[j]==-1){
//                        A=j;
//                    }
//                    else if (dummy[j]==-latticedim*(latticedim-1)){
//                        D=j;
//                    }
//                    else if (dummy[j]==1){
//                        B=j;
//                    }
//                    else if (dummy[j]==-latticedim){
//                        C=j;
//                    }
//                    else {
//                        cout << "FAILS" << endl;
//                    }
//                }
//            }
//            else if (dummysum == -latticedim-1){
//    //            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 7 "  << endl;
//                for (int j=0;j<4;j++){
//                    if (dummy[j]==-1){
//                        A=j;
//                    }
//                    else if (dummy[j]==-latticedim*(latticedim-1)){
//                        D=j;
//                    }
//                    else if (dummy[j]==1-latticedim){
//                        B=j;
//                    }
//                    else if (dummy[j]==-latticedim){
//                        C=j;
//                    }
//                    else {
//                        cout << "FAILS" << endl;
//                    }
//                }
//            }
//        }

    }

    RenderWindow app(VideoMode(65+scale*latticeDim*2.5*390/16, 55+scale*latticeDim*2.5*390/16), "The Pipe Puzzle!");

    cout << endl << endl << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Render"<< endl;
    Texture t4,t5;
    t4.loadFromFile("/home/oli/CLionProjects/buildgame/images/pipes.png");
    t4.setSmooth(true);
    t5.loadFromFile("/home/oli/CLionProjects/buildgame/images/unpipes.png");
    t5.setSmooth(true);

    Sprite sPipe(t4);
    sPipe.setOrigin(27,27);
    Sprite sUnpipe(t5);
    sUnpipe.setOrigin(27,27);

    while (app.isOpen()) {
        Event e{};
        Vector2i selectednode = {0, 11};
//        int ts = 54; //lattice grid size
        int pngts = 54; //lattice grid size

        int ts = pngts*scale;

//        sf::Vector2f offset(65, 55);
        sf::Vector2f offset(65, 55);


        while (app.pollEvent(e)) {


            if (e.type == Event::Closed) {
                converged = true;
                app.close();
            } else converged = false;

            if (e.type == Event::MouseButtonPressed) {
                if (e.key.code == Mouse::Left) {
                    Vector2i pos = Mouse::getPosition(app) + Vector2i(ts / 2, ts / 2) - Vector2i(offset);
                    pos /= ts;
                    if (isOut(pos, latticeDim)) continue;
                    cout << "Registered Left Click at " << pos.x <<"  " << pos.y << endl;
                    cell(pos).orientation++;
                    cell(pos).orientation = cell(pos).orientation % 4;
                    cell(pos).rotate();
                    int original_orientation = cell(pos).orientation;
                    selectednode = pos;

                    //noderef
                    int node_ref = int(N * pos.y + pos.x);
                    cout << "--associated with node : " << node_ref << endl;
                    cout << "--orientation : " << cell(pos).orientation << endl;
                    nodes[node_ref].clockwiseCnxs(latticeDim);

                    pipes[node_ref].setCnxs(nodes[node_ref].cnxs);
                    pipes[node_ref].cnxstodirs();
                    pipes[node_ref].setOrientation();
                    pipes[node_ref].setAngle();

                    int final_orientation = cell(pos).orientation;
                    cout << "Orientations" << endl;
                    cout << original_orientation << "   " << final_orientation << endl;
                    cout << "IDS" << endl;
                    cout << "node : " << nodes[node_ref].id << endl;
                    cout << "pipe : " << pipes[node_ref].id << endl;
                    cout << endl;
//                    cout << "CONNECTIONS : " << endl;
//                    cout << "Node : " << nodes[node_ref].cnxs << endl;
//                    cout << "Pipes: " << cell(pos).cnxs << endl;
                    cout << "--rotated node" << endl;
                    for (int i = 0; i < N*N; i++) {
                        pipe &p = pipes[i];
                        p.on = 0;
                    }

                    pipe &p = pipes[pos.x+N*pos.y];
                    //p.rotate();
                    int cumulativeEnergy = 0;
                    for (int i = 0; i < N * N; i++) {
                        int x = int(i % N);
                        int y = (int((i - x) / N));
                        pipe &p = pipes[i];
                        std::string s = "";
                        if (x == pos.x && y == pos.y) {
                            cout << "String : " << p.orientationstring() << endl;
                            cout << "Orient : " << p.orientation << endl;
                            cout << "Node E : " << p.energynode(pipes, latticeDim) << endl;
                        }
                        cumulativeEnergy += p.energynode(pipes, latticeDim);
//                        for (int j = 0; j < originalsize; j++) {
////                            cout << j << "  :  " << energy << endl;
//                            int x_trial;
//                            int y_trial;
//                            ////
//                            if (p.dirs[j].x < 0) {
//                                x_trial = x0;
//                                y_trial = y;
//
////                                cout << x_trial <<"," << y_trial << endl;
//                                pipe &q = grid[x_trial][y_trial];
//                                bool connect = false;
//                                int newsize = q.dirs.size();
//                                for (int k = 0; k < newsize; k++) {
//                                    if (q.dirs[k] == Right) connect = true;
//                                }
//
//                                if (x == pos.x && y== pos.y) cout << "LEFT : " << connect << endl;
//                                if (connect == false) energy += 1; // not true! need unsated cnxs
//                            } else if (p.dirs[j].x > 0) {
//                                x_trial = x1;
//                                y_trial = y;
////                                cout << x_trial <<"," << y_trial << endl;
//                                pipe &q = grid[x_trial][y_trial];
//                                bool connect = false;
//                                int newsize = q.dirs.size();
//                                for (int k = 0; k < newsize; k++) {
//                                    if (q.dirs[k] == Left) connect = true;
//                                }
//                                if (x == pos.x && y== pos.y) cout << "RIGHT : " << connect << endl;
//                                if (connect == false) energy += 1;
//                            } else {
//                                x_trial = x;
//                                if (p.dirs[j].y < 0) {
//                                    y_trial = y0;
////                                    cout << x_trial <<"," << y_trial << endl;
//                                    pipe &q = grid[x_trial][y_trial];
//                                    bool connect = false;
//                                    int newsize = q.dirs.size();
//                                    for (int k = 0; k < newsize; k++) {
//                                        if (q.dirs[k] == Down) connect = true;
//                                    }
//                                    if (x == pos.x && y== pos.y) cout << "TOP : " << connect << endl;
//                                    if (connect == false) energy += 1;
//                                } else if (p.dirs[j].y > 0) {
//                                    y_trial = y1;
////                                    cout << x_trial <<"," << y_trial << endl;
//                                    pipe &q = grid[x_trial][y_trial];
//                                    bool connect = false;
//                                    int newsize = q.dirs.size();
//                                    for (int k = 0; k < newsize; k++) {
//                                        if (q.dirs[k] == Up) connect = true;
//                                    }
//                                    if (x == pos.x && y== pos.y) cout << "BOTTOM : " << connect << endl;
//                                    if (connect == false) energy += 1;
//                                }
//                            }
//
//                        }
//
                    }
                    cout << "Energy : " << cumulativeEnergy << endl;
                    cout << endl;
                }

                if (e.key.code == Mouse::Right) {
                    Vector2i pos = Mouse::getPosition(app) + Vector2i(ts / 2, ts / 2) - Vector2i(offset);
                    pos /= ts;
                    if (isOut(pos, latticeDim)) continue;

                    cell(pos).bendorstraighten();
                    cout << "sorted pipe" << endl;
                    int node_ref = int(N*pos.y+pos.x);
                    cout << "NODE REF = " << node_ref << endl;
                    if (cell(pos).dirs.n==4){
                        nodes[node_ref].cnxs = nodes[node_ref].allowedCnxs;
                    }
                    else if (cell(pos).dirs.n==3){
                        nodes[node_ref].cnxs.setSize(3);
                        nodes[node_ref].setOrientation(cell(pos).orientation, latticeDim);
                    }

//                    nodes[node_ref].bendorstraight(cell(pos).dirs.n, threecoordcount, fourcoordcount, threecoord, fourcoord);
                    cout << "sorted node" << endl;
                    cell(pos).setCnxs(nodes[node_ref].cnxs);
                    cell(pos).cnxstodirs();
                    cell(pos).setOrientation();
                    cell(pos).setAngle();

//                    cout << "CONNECTIONS : " << endl;
//                    cout << "Node : " << nodes[node_ref].cnxs << endl;
//                    cout << "Pipes: " << cell(pos).cnxs << endl;

                    int newfourcoordcount = 0;
                    int newthreecoordcount = 0;
                    int newtwocoordcount = 0;
                    int dummy;
                    for (int i = 0; i < N * N; i++) {
                        int x = i % N;
                        int y = int((i - x) / N);
                        dummy = int(cell({x, y}).dirs.n);
                        if (dummy == 2) {
                            newtwocoordcount += 1;
                        } else if (dummy == 4) {
                            newfourcoordcount += 1;
                        } else if (dummy == 30 or dummy == 31 or dummy == 3) {
                            newthreecoordcount += 1;
                        } else {
                            cout << "UNKNOWN VALUE" << endl;
                        }
                    }
                    cout << "               N+1 :  " << newfourcoordcount << endl;
                    cout << "               N   :  " << newthreecoordcount << endl;
                    cout << "               N-2 :  " << newtwocoordcount << endl;
                    cout << "               --------------------------" << endl;
                    cout << "               TOT :  " << newfourcoordcount + newthreecoordcount + newtwocoordcount
                         << endl;


                    cout << "ENERGY  :  " << pipeenergytot(latticeDim) << endl;

                }
            }

            if (e.type == Event::KeyPressed) {
                if (e.key.code == sf::Keyboard::Tab) {
//                    ofstream cnxsFile("gameoutput_"+to_string(energytot())+".dat");
//                    if (cnxsFile.is_open()){
//                        for(int i=0; i<N*N; ++i) {
//                            int x = i%N;
//                            int y = int((i-x)/N);
//                            for (int j=0;j<cell({x,y}).dirs.size();j++){
//                                int x0 = x+cell({x,y}).dirs[j].x;
//                                int y0 = y+cell({x,y}).dirs[j].y;
//                                //cout <<"    ("<<x0<<","<<y0<<")" << endl;
//                                if (x0<0) x0+=(N);
//                                else if (x0>N-1) x0-=(N);
//                                if (y0<0) y0+=(N);
//                                else if (y0>N-1) y0-=(N);
//                                int ref = x0 + N*y0;
//                                cout << "    "<<i<<"    "<<ref<<endl;
//                                cnxsFile<<to_string(i)+"    "+to_string(ref)+"\n";
//                            }
//                        }
//
//                    }

                    ofstream cnxsFile("gameoutput.dat");
                    if (cnxsFile.is_open()) {
                        cnxsFile << to_string(N * N);
                        for (int i = 0; i < N * N; ++i) {
                            int x = i % N;
                            int y = int((i - x) / N);
                            cnxsFile << "\n";
                            for (int j = 0; j < cell({x, y}).dirs.n; j++) {
                                int x0 = x + cell({x, y}).dirs[j].x;
                                int y0 = y + cell({x, y}).dirs[j].y;
                                //cout <<"    ("<<x0<<","<<y0<<")" << endl;
                                if (x0 < 0) x0 += (N);
                                else if (x0 > N - 1) x0 -= (N);
                                if (y0 < 0) y0 += (N);
                                else if (y0 > N - 1) y0 -= (N);
                                int ref = x0 + N * y0;
                                //cout << "    "<<i<<"    "<<ref<<endl;
                                cnxsFile << setw(20) << left << to_string(ref) + "    ";
                            }
                        }


                    }

                }
            }

            if (e.type == Event::KeyPressed) {
                if (e.key.code == sf::Keyboard::Enter) {
                    //begin an mc round of 1000 steps
                    float temp;
                    int step = 0;
                    int initialenergy = pipeenergytot(latticeDim);
                    while (pipeenergytot(latticeDim) > 0.5 * initialenergy && step < 10000) {
//                    while (energytot()>0){

                        int mcnode = rand() % (N * N);
                        int numberrotations = rand() % 3;

                        int x = mcnode % N;
                        int y = int((mcnode - x) / N);
                        int node_ref = int(N * y + x);
                        Vector2i neighbour1 = putin({x + 1, y}, latticeDim);
                        Vector2i neighbour2 = putin({x - 1, y}, latticeDim);
                        Vector2i neighbour3 = putin({x, y + 1}, latticeDim);
                        Vector2i neighbour4 = putin({x, y - 1}, latticeDim);
                        bool run = true;

                        if (pipeenergyVector({x, y}, latticeDim) != 0) {
                            run = true;
                        } else if (pipeenergyVector(neighbour1, latticeDim) != 0 or pipeenergyVector(neighbour2, latticeDim) != 0 or
                                   pipeenergyVector(neighbour3, latticeDim) != 0 or pipeenergyVector(neighbour4, latticeDim) != 0) {
                            run = true;
                        }

//                        cout << x << "," << y << endl;
                        if (run) {

                            if (step % 100 == 0)
                                cout << "     STEP : " << step << "      " << pipeenergytot(latticeDim) << " / " << initialenergy
                                     << endl;
                            step++;
                            Vector2i pos = {x, y};
                            int e0 = pipeenergytot(latticeDim);

                            for (int j = 0; j < numberrotations + 1; j++) {
                                cell(pos).orientation++;
                                cell(pos).orientation = cell(pos).orientation % 4;
                                cell(pos).rotate();
                                nodes[node_ref].clockwiseCnxs(latticeDim);
                            }

                            int e1 = pipeenergytot(latticeDim);

                            if (e1 > e0) {
                                float rando = ((double) rand() / (RAND_MAX));
                                if (e0 < 6 && e0 > 4) temp = 0.001;
                                else if (e0 <= 4) temp = 0.0001;
                                else temp = 0.01;
//                                cout << "RANDO : " << rando << "  vs : " << exp((e0-e1)/temp)<< endl;
                                if (temp == 0) {
                                    for (int j = 0; j < (3 - numberrotations); j++) {
                                        cell(pos).orientation++;
                                        cell(pos).rotate();
                                        nodes[node_ref].clockwiseCnxs(latticeDim);
                                    }
                                } else if (rando > exp((e0 - e1) / temp)) {
                                    for (int j = 0; j < (3 - numberrotations); j++) {
                                        cell(pos).orientation++;
                                        cell(pos).rotate();
                                        nodes[node_ref].clockwiseCnxs(latticeDim);
                                    }

                                }
                            }
                        }
                    }

                    int newfourcoordcount = 0;
                    int newthreecoordcount = 0;
                    int newtwocoordcount = 0;
                    int dummy;
                    for (int i = 0; i < N * N; i++) {
                        int x = i % N;
                        int y = int((i - x) / N);
                        dummy = int(cell({x, y}).dirs.n);
                        if (dummy == 2) {
                            newtwocoordcount += 1;
                        } else if (dummy == 4) {
                            newfourcoordcount += 1;
                        } else if (dummy == 30 or dummy == 31 or dummy == 3) {
                            newthreecoordcount += 1;
                        } else {
                            cout << "UNKNOWN VALUE" << endl;
                        }
                    }
                    cout << "               N+1 :  " << newfourcoordcount << endl;
                    cout << "               N   :  " << newthreecoordcount << endl;
                    cout << "               N-2 :  " << newtwocoordcount << endl;
                    cout << "               --------------------------" << endl;
                    cout << "               TOT :  " << newfourcoordcount + newthreecoordcount + newtwocoordcount
                         << endl;


                    cout << "ENERGY  :  " << pipeenergytot(latticeDim) << endl;
                }
            }
        }
//        cout << "reached app clear" << endl;
        app.clear();
        //app.draw(sBackground);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//        cout << "-----> #sanity check ... ";
        for (int i=0;i<N*N;++i){
            for (int j=0;j<nodes[i].cnxs.n;++j){
                if (nodes[i].cnxs[j] != pipes[i].cnxs[j]){
                    cout << "for fucks sake" << endl;
                    cout << nodes[i].cnxs[j] << "  " << pipes[i].cnxs[j] << endl;
                    cout << endl;
                }
            }
            if (nodes[i].cnxs.n!=pipes[i].cnxs.n) cout << nodes[i].cnxs.n << "  " << pipes[i].cnxs.n << endl;
        }
//        cout << "completed!" << endl;

//        for (int j=0;j<display.size();j++){
//            int i=display[j];
//            cout << j << ":" << i << endl;
        int x;
        int y;
        for (int i = 0; i < N * N; i++) {
            x = ((i) % N);
            y = (int((i - x) / N));
            //0,12
//        for (int i=0;i<2;i++){
//            if (i==1){
//                x = selectednode.x;
//                y = selectednode.y;
//            }
//            else {
//                pipe &q = grid[selectednode.x][selectednode.y];
//                x = selectednode.x + q.dirs[i].x;
//                y = selectednode.y + q.dirs[i].y;
//                if (x<0) x+=N;
//                else if (x>N) x-=N;
//                if (y<0) y+=N;
//                else if (y>N) y-=N;
//            }
            pipe &p = pipes[i];

//            p.angle += 5;


            p.angle += 5;
            if (p.angle > p.orientation * 90) p.angle = p.orientation * 90;

            int dummy = p.dirs.n;

            std::string s = p.orientationstring();

            if (p.energynode(pipes, latticeDim) == 0) {
                //                cout << "line 523" << endl;

                if (dummy == 2) {
                    if (s == "1010" or s == "0101") {
                        sPipe.setTextureRect(IntRect(0, 0, pngts, pngts));
                    }
                    else {
                        sPipe.setTextureRect(IntRect(2 * pngts, 0, pngts, pngts));
                    }
                }
                else if (dummy == 3) {
                    sPipe.setTextureRect(IntRect(3 * pngts, 0, pngts, pngts));
                }
                else if (dummy == 4) {
                    sPipe.setTextureRect(IntRect(pngts, 0, pngts, pngts));

                }
                //                cout << "line 531" << endl;
                sPipe.setScale(scale, scale);
                sPipe.setRotation(p.angle);
                sPipe.setPosition(x * ts, y * ts);
                sPipe.move(offset);
                //                cout << "line 536" << endl;
                app.draw(sPipe);
            } else {
                if (dummy == 2) {
                    if (s == "1010" or s == "0101") {
                        sUnpipe.setTextureRect(IntRect(0, 0, pngts, pngts));
                    }
                    else {
                        sUnpipe.setTextureRect(IntRect(2 * pngts, 0, pngts, pngts));
                    }
                }
                else if (dummy == 3) {
                    sUnpipe.setTextureRect(IntRect(3 * pngts, 0, pngts, pngts));
                }
                else if (dummy == 4) {
                    sUnpipe.setTextureRect(IntRect(pngts, 0, pngts, pngts));
                }

                sUnpipe.setScale(scale, scale);
                sUnpipe.setRotation(p.angle);
                sUnpipe.setPosition(x * ts, y * ts);
                sUnpipe.move(offset);
                app.draw(sUnpipe);
            }


        }
        app.display();

    }



    //Set up status
    VecF<int> opt(2);
    opt[0]=converged;
    opt[1]=10;

    if(converged){
        if(calcRings) {
            int ringStatus = findRings();
            if (ringStatus == 1) opt[0] = 2;
            else if (ringStatus == 2) opt[0] = 3;
            if (opt[0] == 1) if (calcEnvs) findEnvironments();
        }
        if(calcChains){
            findChains();
        }
    }
//    writeNetwork("./output/mc",iterations);

    return opt;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        for(int i=0;i<N;i++) {
//            for (int j = 0; j < N; j++) {
//                pipe &p = grid[j][i];
//                int kind = p.dirs.size();
////                p.dirs = {};
//                cout << i<<","<<j<<"    "<<p.coordinaitionno() - kind << endl;
//                if (kind == 2 && p.dirs[0] == -p.dirs[1]) kind = 0;
//                if (kind == 1) kind = 3;
//                p.angle += 5;
//                if (p.angle > p.orientation * 90) p.angle = p.orientation * 90;
//
//                sPipe.setTextureRect(IntRect(ts * kind, 0, ts, ts));
//                sPipe.setRotation(p.angle);
//                sPipe.setPosition(j * ts, i * ts);
//                sPipe.move(offset);
//                app.draw(sPipe);
//
////                 if (kind==1)
////                     { if (p.on) sComp.setTextureRect(IntRect(53,0,36,36));
////                       else sComp.setTextureRect(IntRect(0,0,36,36));
////                       sComp.setPosition(j*ts,i*ts);sComp.move(offset);
////                       app.draw(sComp);
////                     }
//            }
//        }
//
////        app.draw(sServer);
//        app.display();
}

VecF<int> Lattice::generate(int maxIterations, double temperature, int mcsw) {

    //Monte Carlo search for fully coordinated lattice
    mcswitch = mcsw;

    int orderedfourcoord = fourcoordcount*fractionorderedH;


    generater2E();

    cout << "ORIGINAL ENERGY  : "<< energy << endl;
    cout << "ORIGINAL R2      : "<< r2 << endl;
    //energy -= 0;


    //Monte Carlo optimisation
    int iterations = 0;
    int cutoffIterations = 0;
    int cutoff = nodes.n*100;

    converged = false;

    //following line applies speedup if n+1 coordinate sites do not change on rotation (ie maximised)
//    uniform_int_distribution<int> randNode(orderedfourcoord,nodes.n-1);
//    uniform_int_distribution<int> randNode(0,threecoordcount-1);
    uniform_int_distribution<int> randNode(0,nodes.n-1);
    uniform_real_distribution<float> mcrand(0.0,1.0);

    int energyE;

    for(;;) {
    //following line applies speedup if n+1 coordinate sites do not change on rotation (ie maximised)
        int id = randNode(mtGen);
        int idcoord = randNode(mtGen);
        int coordoriginal;
        int coord;
        coord = ncoordsoriginal[id];

        ////////////////////////////original cnxs////////////////////////////////////////////////////////////////
        int currCnd=nodes[id].getReciprocalCnxs(nodes).n;
        VecR<int> currCnxs=nodes[id].cnxs;

        ////////////////////////////randomises cnxs for trial////////////////////////////////////////////////////
        if (coord != 4) {
            if (!pattern) {
                nodes[id].randomCnxs(coord, mtGen, latticeDim);
            }
            else {
                if (randRL(mtGen) == 0) {
                    nodes[id].randomCnxs(coord, patternR, mtGen);
                }
                else {
                    nodes[id].randomCnxs(coord, patternL, mtGen);
                }
            }
        }
        ////////////////////////////randomising id can only affect neightbouring e///////////////////////////////

        int trialCnd=nodes[id].getReciprocalCnxs(nodes).n;
        int deltaCnd=trialCnd-currCnd;

//        if (trialCnd > currCnd){
//           cout << currCnd << "-->" << trialCnd << "/" << nodes[id].cnxs.n << "    " << energy << "    " << energy - 2*deltaCnd << "    " << endl;
//        }
        if(deltaCnd<0){
            if(rand01(mtGen)<exp(deltaCnd/temperature)){
                energy-=2*deltaCnd;
                cutoffIterations=0;
            }
            else{
                nodes[id].setCnxs(currCnxs);
                cutoffIterations+=1;
            }
        }
        else {
            energy -= 2 * deltaCnd;
            cutoffIterations = 0;
        }


        if(energy==0){
            checkCoords();
            if (check4coord == fourcoordcount && check3coord == threecoordcount && check2coord == twocoordcount ) {
                converged=true;
                cout << "CONVERGED" <<endl;
                break;
            }
            else{
                cutoffIterations = 0;
                cout << "Coordinate ratios not consistent" << endl;
            }
        }
        else if(cutoffIterations>=cutoff) {
            cout << "REACHED CUTOFF ITERATIONS" << endl;
//            cout << "Trying md ..." << endl;
//            cutoffIterations=0;
//            iterations=0;
//            for(;;){
//                id = selectgoodnode();
//                generater2E();
//                int energy0 = energy;
//                int r20 = r2;
//                optimumCnxs(id);
//                generater2E();
//                int energy1 = energy;
//                int r21= r2;
//                if (energy1>=energy0){
//                    cutoffIterations+=1;
//                }
//                if (energy0!=energy1){
//                    cout << iterations << "  " << energy0 << " --> " << energy1 << endl;
//                }
//                if (energy==0){
//                    converged=true;
//                    cout << "CONVERGED" <<endl;
//                    break;
//                }
//                else if(cutoffIterations>=cutoff){
//                    cout << "REACHED CUTOFF ITERATIONS" << endl;
//                    break;
//                }
//                else if(iterations>=maxIterations){
//                    cout << "REACHED MAX ITERATIONS" << endl;
//                    break;
//                }
//                ++iterations;
//            }

            break;
        }
        else if(iterations>=maxIterations) {
            cout << "REACHED MAX ITERATIONS" << endl;
//            cout << "Trying md ..." << endl;
//            cutoffIterations=0;
//            iterations=0;
//            generater2E();
//            for(;;) {
//                id = selectgoodnode();
//                generater2E();
//                int energy0 = energy;
//                int r20 = r2;
//                optimumCnxs(id);
//                generater2E();
//                int energy1 = energy;
//                int r21= r2;
//                if (energy1>=energy0){
//                    cutoffIterations+=1;
//                }
//                if (energy0!=energy1){
//                    cout << iterations << "  " << energy0 << " --> " << energy1 << endl;
//                }
//                if (energy==0){
//                    converged=true;
//                    cout << "CONVERGED" <<endl;
//                    break;
//                }
//                else if(cutoffIterations>=cutoff){
//                    cout << "REACHED CUTOFF ITERATIONS" << endl;
//                    break;
//                }
//                else if(iterations>=maxIterations){
//                    cout << "REACHED MAX ITERATIONS" << endl;
//                    break;
//                }
//                ++iterations;
//                cout << "ID : " << id << endl;
//                generater2E();
//                int energy0 = energy;
//                int r20 = r2;
//                cout << "e0 : " << energy0 << "    r0 : " << r20 << endl;
//                optimumCnxs(id);
//                generater2E();
//                int energy1 = energy;
//                int r21= r2;
//                cout << "e1 : " << energy1 << "    r1 : " << r21 << endl;
//                if (energy1>=energy0){
//                    cutoffIterations+=1;
//                }
//                if (energy1>energy0){
//                    cout << "ENERGY ERROR" << endl;
//                }
//                if (energy0!=energy1){
//                    cout << iterations << "  " << energy0 << " --> " << energy1 << endl;
//                }
//                if (energy==0){
//                    converged=true;
//                    cout << "CONVERGED" <<endl;
//                    break;
//                }
//                else if(cutoffIterations>=cutoff){
//                    cout << "REACHED CUTOFF ITERATIONS" << endl;
//                    break;
//                }
//                else if(iterations>=maxIterations){
//                    cout << "REACHED MAX ITERATIONS" << endl;
//                    break;
//                }
//                ++iterations;
//
//            }

            break;
        }
        else if(energy<-100){
            cout << "ENERGY PASSED ZERO ERR" << endl;
            break;
        }

        ++iterations;
//        cout << "CUTOFF  :  " << cutoff << endl;
//        cout << "MAX ITERATIONS  :  " << maxIterations << endl;
//        if(iterations%1==0) writeNetwork("./vis",iterations);
    }

    //Set up status
    VecF<int> opt(2);
    opt[0]=converged;
    opt[1]=iterations;

    if(converged){
        if(calcRings) {
            int ringStatus = findRings();
            if (ringStatus == 1) opt[0] = 2;
            else if (ringStatus == 2) opt[0] = 3;
            if (opt[0] == 1) if (calcEnvs) findEnvironments();
        }
        if(calcChains){
            findChains();
        }
    }
//    writeNetwork("./output/mc",iterations);

    return opt;
}


int Lattice::generate(VecR<int> &pairA, VecR<int> &pairB) {
    //Generate from defined list of edges to break
    cout << "WELL WE're using the PAIR GENERATE THEN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    //Loop over edges and remove connections
    for(int i=0; i<pairA.n; ++i){
        int idA=pairA[i];
        int idB=pairB[i];
        nodes[idA].breakCnx(idB);
        nodes[idB].breakCnx(idA);
    }

    //Find rings
    findRings();

    //Find environments
    if(calcEnvs) findEnvironments();

    //Calculate energy
    int energy=nodes.n*nodeCnd;
    for(int i=0; i<nodes.n; ++i){
        energy -= nodes[i].getReciprocalCnxs(nodes).n;
    }

    return energy;
}


int Lattice::findRings() {

    //Merge rings
    for(int i=0; i<nodes.n; ++i){
        //Get broken connections
        VecR<int> vacantCnxs=nodes[i].getVacantCnxs();
//        cout << "VACANT CNXS : " << vacantCnxs << endl;
        for(int j=0; j<vacantCnxs.n; ++j){
            if(i<vacantCnxs[j]){
                 //prevent double counting
                 //Get ids of rings that nodes share
                 VecR<int> shared=nodes[i].getSharedRings(nodes[vacantCnxs[j]]);
                 VecR<int> rIds(0,shared.n);
                 for(int k=0; k<shared.n; ++k){
                     if(rings[shared[k]].checkEdge(i,vacantCnxs[j])) rIds.addValue(shared[k]);
                 }
                 if(rIds.n!=2){
//                     writeCrds("./split");
//                     writeNetwork("./split",i);
                     return 1; //sample split in two so fails
                 }
                 int rid0=rIds[0];
                 int rid1=rIds[1];
//                 cout << "RINGS TO MERGE : " << rid0 << "  " << rid1 << endl;
                 //Generate new ring by merging shared rings
                 Ring ring=rings[rid0].merge(rings[rid1],i,vacantCnxs[j]);
                 int rId2=ring.id;
                 //Deactivate old rings and change ids of associated nodes
                 rings[rid0].clear();
                 rings[rid1].clear();
                 for(int k=0; k<ring.nodes.n; ++k){
                     nodes[ring.nodes[k]].swapRing(rid0,rId2);
                     nodes[ring.nodes[k]].swapRing(rid1,rId2);
                 }
                 rings.addValue(ring);
            }
        }
    }

    //Find ring adjacencies
    int nid0,nid1;
    for(int i=0; i<nodes.n; ++i){
        nid0=i;
        for(int j=0; j<nodes[i].cnxs.n; ++j){
            nid1=nodes[i].cnxs[j];
            if(nid0<nid1){
                VecR<int> shared=nodes[nid0].getSharedRings(nodes[nid1]);
                VecR<int> rIds(0,shared.n);
                for(int k=0; k<shared.n; ++k){
                    if(rings[shared[k]].checkEdge(nid0,nid1)) rIds.addValue(shared[k]);
                }
                if(rIds.n==1){//self interaction
                    rings[rIds[0]].addCnx(rIds[0]);
                    rings[rIds[0]].addCnx(rIds[0]);
                }
                else if(rIds.n==2){
                    rings[rIds[0]].addCnx(rIds[1]);
                    rings[rIds[1]].addCnx(rIds[0]);
                }
                else{//should not occur
                    return 2;
                }
            }
        }
    }

    //Check number of adjacencies match number of nodes
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            if(rings[i].cnxs.n!=rings[i].nodes.n) return 2;
        }
    }

    //Calculate coordinates
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            VecF<double> x(rings[i].nodes.n),y(rings[i].nodes.n);
            for(int j=0; j<rings[i].nodes.n; ++j){
                x[j] = nodes[rings[i].nodes[j]].crd[0];
                y[j] = nodes[rings[i].nodes[j]].crd[1];
            }
            VecF<double> com=periodicCentreOfMass(x,y);
            //rings[i].setCrd(com);
        }
    }

    return 0;
}


int Lattice::findChains() {
    //Find chains of 2-coordinate nodes

    //Initialise chain vector
    chains=VecR<Chain>(0,nodes.n/2+1);

    //Find chains by following node connections until reach original node
    VecF<bool> nodeInc(nodes.n);
    nodeInc=false;
    for(int i=0; i<nodes.n; ++i){
        if(!nodeInc[i]){
            Chain chain(nodes.n);
            chain.addNode(i);
            int nId=i;
            int nIdPrev=nodes[nId].cnxs[0];
            for(;;){
                int nIdNext;
                if(nodes[nId].cnxs[0]==nIdPrev) nIdNext=nodes[nId].cnxs[1];
                else nIdNext=nodes[nId].cnxs[0];
                if(nIdNext==i) break;
                else{
                    chain.addNode(nIdNext);
                    nIdPrev=nId;
                    nId=nIdNext;
                }
            }
            for(int j=0; j<chain.nodes.n; ++j) nodeInc[chain.nodes[j]]=true;
            chains.addValue(chain);
        }
    }

    return 0;
}


void Lattice::chainAnalysis(VecR<int> &chainLengths) {
    //Analyse chains

    //Get chain lengths
    for(int i=0; i<chains.n; ++i){
        chainLengths.addValue(chains[i].nodes.n);
    }
}


VecF<double> Lattice::periodicCentreOfMass(VecF<double> &x, VecF<double> &y) {
    //Find centre of mass accounting for periodic boundary conditions

    VecF<double> thetaX,thetaY;
    thetaX = x*(2*M_PI/periodicity[0]);
    thetaY = y*(2*M_PI/periodicity[1]);
    VecF<double> alphaX,betaX,alphaY,betaY;
    alphaX = vCos(thetaX);
    betaX = vSin(thetaX);
    alphaY = vCos(thetaY);
    betaY = vSin(thetaY);
    double avAlphaX=vMean(alphaX);
    double avBetaX=vMean(betaX);
    double avAlphaY=vMean(alphaY);
    double avBetaY=vMean(betaY);
    double phiX = atan2(-avBetaX,-avAlphaX)+M_PI;
    double phiY = atan2(-avBetaY,-avAlphaY)+M_PI;
    VecF<double> com(2);
    com[0]=phiX*periodicity[0]/(2*M_PI);
    com[1]=phiY*periodicity[1]/(2*M_PI);

    return com;
}


void Lattice::findEnvironments() {
    //Assign environment type to each original ring in lattice


    if(latticeCode=="sq3"){
        //Square 3-coordinate lattice
        int dim=sqrt(nodes.n);
        //Find environments
        for(int i=0; i<nodes.n; ++i){
            //Active means all sides complete
            if(rings[i].active) rings[i].environment=0;
            else{
                //Otherwise determine missing sides
                VecF<int> ids(4);
                int row=floor(i/dim);
                ids[0]=row*dim+(i+dim-1)%dim;
                ids[1]=((row+1)*dim+(i+dim-1)%dim)%nodes.n;
                ids[2]=((row+1)*dim+i%dim)%nodes.n;
                ids[3]=row*dim+i%dim;
                VecF<int> vacancies(4);
                for(int j=0,k=1; j<4; ++j, ++k){
                    int id0=ids[j];
                    int id1=ids[k%4];
                    vacancies[j]=!vContains(nodes[id0].cnxs,id1);
                }
                //Single vacancy
                int nVac=vSum(vacancies);
                if(nVac==1){
                    for(int j=0; j<4; ++j){
                        if(vacancies[j]==1) rings[i].environment=j+1;
                    }
                }
                else if(nVac==2){
                    if(vacancies[0]==1 && vacancies[2]==1) rings[i].environment=5;
                    if(vacancies[1]==1 && vacancies[3]==1) rings[i].environment=6;
                }
            }
        }
    }
}

VecF<double> Lattice::networkAnalysis(VecF<int> &k, VecF<int> &pk, VecF< VecF<int> > &ejk) {
    //Calculate ring statistics and assortativity

    //Calculate ring statistics and adjacency matrix
    pk=VecF<int>(k.n);
    ejk=VecF< VecF<int> >(k.n);
    for(int i=0; i<k.n; ++i) ejk[i]=VecF<int>(k.n);
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            int ki=rings[i].getNumCnxs();
            for(int j=0; j<rings[i].cnxs.n; ++j){
                int kj=rings[rings[i].cnxs[j]].getNumCnxs();
                ++ejk[ki][kj];
            }
            ++pk[ki];
        }
    }

    VecF<double> res=calculateNetworkProperties(k,pk,ejk);
    return res;
}


VecF<double> Lattice::calculateNetworkProperties(VecF<int> &k, VecF<int> &pk, VecF<VecF<int> > &ejk) {
    //Calculate ring statistics, assortativity and Aboav-Weaire parameters

    //Normalisation
    double pNorm=vSum(pk),eNorm=0.0;
    for(int i=0; i<k.n; ++i) eNorm+=vSum(ejk[i]);

    //Moments and variance
    double k1,k2,k3;
    k1=vSum(k*pk)/pNorm;
    k2=vSum(k*k*pk)/pNorm;
    k3=vSum(k*k*k*pk)/pNorm;
    double mu_2 = k2-k1*k1;

    //Assortativity
    double assort=0.0;
    for(int i=0; i<k.n; ++i){
        for(int j=0; j<k.n; ++j){
            assort+=i*j*ejk[i][j];
//            if(i>=3 && i<=12){
//                if(j>=3 && j<=12) cout<<i<<" "<<j<<" "<<ejk[i][j]<<endl;
//            }
        }
    }
    assort/=eNorm;
    assort=k1*k1*assort-k2*k2;
    assort/=k1*k3-k2*k2;

    //Aboav-Weaire
    VecF<double> qk(k.n),mk(k.n);
    for(int i=0; i<k.n; ++i){
        qk[i]=vSum(ejk[i]);
        mk[i]=vSum(k*ejk[i]);
    }
    VecR<double> awX(0,k.n),awY(0,k.n);
    for(int i=0; i<k.n; ++i){
        if(qk[i]>0){
            awX.addValue(k1*(k[i]-k1));
            awY.addValue(k[i]*mk[i]/qk[i]);
        }
    }
    VecR<double> aw=vLinearRegression(awX,awY);

    //Get principle ring size
    int kMain;
    if(nodeCnd==3) kMain=6;
    else if(nodeCnd==4) kMain=4;
    else if(nodeCnd==5) kMain=3;

    //Return summary
    VecF<double> summary(7);
    summary[0]=pk[kMain]/pNorm;
    summary[1]=k1;
    summary[2]=mu_2;
    summary[3]=assort;
    summary[4]=1.0-aw[0];
    summary[5]=aw[1]-k1*k1;
    summary[6]=aw[2];

    return summary;
}


VecF<int> Lattice::getEnvironments() {
    //Get ring environments

    //Get environment for each ring
    VecF<int> envs(nodes.n);
    for (int i = 0; i < nodes.n; ++i){
        envs[i] = rings[i].environment;
    }

    return envs;
}


string Lattice::getEnvironmentCode() {
    //Get ring environment code

    //Get environment for each ring
    string envCode="";
    for (int i = 0; i < nodes.n; ++i){
        envCode += to_string(rings[i].environment);
    }

    return envCode;
}


void Lattice::rdfAnalysis(VecF<double> &latticeRDF, VecF<double> &dualRDF, double delta) {
    //Calculate rdfs for lattice points and dual

    //Lattice coordinates
    VecF<double> xLat(nodes.n),yLat(nodes.n);
    for(int i=0; i<nodes.n; ++i){
        xLat[i]=nodes[i].crd[0];
        yLat[i]=nodes[i].crd[1];
    }
    calculateRDF(xLat,yLat,latticeRDF,delta);

    //Dual coordinates
    VecF<double> xDual(Ring::totalActive),yDual(Ring::totalActive);
    int j=0;
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            xDual[j]=rings[i].crd[0];
            yDual[j]=rings[i].crd[1];
            ++j;
        }
    }
    calculateRDF(xDual,yDual,dualRDF,delta);
}


void Lattice::calculateRDF(VecF<double> &x, VecF<double> &y, VecF<double> &rdf, double delta) {
    //Calculate rdf

    //Get periodic information
    VecF<double> pbc=periodicity,rpbc(2);
    VecF<double> mic=periodicity/2.0;
    rpbc[0]=1.0/periodicity[0];
    rpbc[1]=1.0/periodicity[1];
    int maxBin;
    double maxDist,maxDistSq;
    if(pbc[0]<=pbc[1]) maxDist=pbc[0]/2.0;
    else maxDist=pbc[1]/2.0;
    maxDistSq=maxDist*maxDist;
    maxBin=floor(maxDist/delta)+1;

    //Calculate pairwise distances and bin
    double xI,yI,b;
    double dx,dy,dSq,d;
    for (int i=0; i<x.n-1; ++i) {
        xI = x[i];
        yI = y[i];
        for (int j=i+1; j<x.n; ++j) {
            dx = xI - x[j];
            dy = yI - y[j];
            dx -= pbc[0] * nearbyint(dx * rpbc[0]);
            dy -= pbc[1] * nearbyint(dy * rpbc[1]);
            dSq = dx * dx + dy * dy;
            if (dSq < maxDistSq) {
                d = sqrt(dSq);
                b = floor(d / delta);
                rdf[b] += 2;
            }
        }
    }
}


void Lattice::skAnalysis(VecF<double> &sk, VecF<int> &multiplicity, double delta, int nMax) {
    //Calculate structure factors for lattice points and dual

    //Lattice coordinates
    VecF<double> xLat(nodes.n),yLat(nodes.n);
    for(int i=0; i<nodes.n; ++i){
        xLat[i]=nodes[i].crd[0];
        yLat[i]=nodes[i].crd[1];
    }
    VecF<double> latticeSk;
//    calculateSk(xLat,yLat,sk,multiplicity,delta,nMax);

    //Dual coordinates
    VecF<double> xDual(Ring::totalActive),yDual(Ring::totalActive);
    int j=0;
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            xDual[j]=rings[i].crd[0];
            yDual[j]=rings[i].crd[1];
            ++j;
        }
    }
    calculateSk(xDual,yDual,sk,multiplicity,delta,nMax);

}


void Lattice::calculateSk(VecF<double> &x, VecF<double> &y, VecF<double> &sk, VecF<int> &multiplicity, double delta, int nMax) {
    //Calculate structure factor

    //Variables
    double kx=2*M_PI/periodicity[0];
    double ky=2*M_PI/periodicity[1];

    //Set up real and imaginary components
    int numN=nMax*2+1;
    VecF< VecF<double> > real(numN),imag(numN),skk(numN);
    for(int nx=-nMax,i=0; nx<=nMax; ++nx,++i){
        real[i]=VecF<double>(numN);
        imag[i]=VecF<double>(numN);
        skk[i]=VecF<double>(numN);
    }
    for(int nx=-nMax,i=0; nx<=nMax; ++nx,++i){
        for(int ny=-nMax,j=0; ny<=nMax; ++ny,++j) {
            VecF<double> v=x*nx*kx+y*ny*ky;
            real[i][j]=vSum(vCos(v));
            imag[i][j]=vSum(vSin(v));
        }
    }

    //Calculate structure factor
    for(int nx=-nMax,i=0; nx<=nMax; ++nx,++i){
        for(int ny=-nMax,j=0; ny<=nMax; ++ny,++j) {
            skk[i][j]=(real[i][j]*real[i][j]+imag[i][j]*imag[i][j])/x.n;
        }
    }

    //Normalise
    double kSq,rDelta=1.0/delta;
    int b;
    for(int nx=-nMax,i=0; nx<=nMax; ++nx,++i){
        for(int ny=-nMax,j=0; ny<=nMax; ++ny,++j) {
            kSq=kx*nx*kx*nx+ky*ny*ky*ny;
            cout<<nx<<" "<<ny<<" "<<skk[i][j]<<endl;
            if(kSq>0) {
                b=floor(sqrt(kSq)*rDelta);
                sk[b]+=skk[i][j];
                ++multiplicity[b];
            }
        }
    }
}


int Lattice::getNumNodes() {
    //Get number of nodes

    return nodes.n;
}


int Lattice::getNumActiveRings() {
    //Get number of active rings

    return Ring::totalActive;
}


VecF<double> Lattice::getPeriodicity() {
    //Get periodicity information

    return periodicity;
}


void Lattice::getEdgeCombinations(VecR<int> &pairA, VecR<int> &pairB, int &n, int &r) {
    //Get unique pairs of edges

    pairA=VecR<int>(0,nodes.n*nodeCnd);
    pairB=VecR<int>(0,nodes.n*nodeCnd);
    int idA,idB;
    for(int i=0; i<nodes.n; ++i){
        idA=i;
        for(int j=0; j<nodes[i].cnxs.n; ++j){
            idB=nodes[i].cnxs[j];
            if(idA<idB){
                pairA.addValue(idA);
                pairB.addValue(idB);
            };
        }
    }
//    cout << pairA << "    " << pairB << endl;
    n = pairA.n;
    r = nodes.n*(latticeCnd-nodeCnd)/2;
}


void Lattice::writeCrds(string prefix) {
    //Write periodicity, node and ring coordinates
    cout << "####################### USING WRITE COORDINATES ###########################" << endl;

    int dummyfour = int(nodes.n * fractionH);
    int dummytwo = int(nodes.n * fractionL);
    int dummythree = nodes.n - dummyfour - dummytwo;

    if (dummytwo % 2 != 0) {
        dummytwo += 1;
        dummythree -= 1;
    }
    if (dummyfour% 2 != 0) {
        dummyfour+= 1;
        dummythree-= 1;
    }
    if ((nodes.n - dummyfour- dummythree- dummytwo) == -1) {
        dummythree-= 1;
    }
    else if ((nodes.n - dummyfour- dummythree- dummytwo) == +1) {
        dummythree += 1;
    }

    int dummylinear;
    if (useLlinear == 1) dummylinear = int(dummytwo*fractionLlinear);
    else dummylinear=dummytwo*0.5;

    string newoutputPrefixfolder;
    string newoutputPrefixfile;
    if (usearray == 1){
        newoutputPrefixfolder = string("./output_"+outputfolder);
        newoutputPrefixfile = string(secondaryinputFile+"_out");
    }
    else if(outStyle==1) {
        newoutputPrefixfolder = "./output/"+outputPrefix +"_4_" + to_string(dummyfour)+"_3_"+to_string(dummythree)+"_2_"+to_string(dummytwo)+"_linear_"+to_string(dummylinear);
        newoutputPrefixfile = "4_" + to_string(dummyfour)+"_3_"+to_string(dummythree)+"_2_"+to_string(dummytwo) +"_linear_"+to_string(dummylinear);
    }
    else {
        newoutputPrefixfolder = "./output";
        newoutputPrefixfile = outputPrefix;
    }

    const char *checknewoutputPrefixfolder = newoutputPrefixfolder.c_str();

    struct stat info;
    if (stat(checknewoutputPrefixfolder, &info) !=0) {
        cout << "LATTICE : CANNOT ACCESS DIRECTORY " << checknewoutputPrefixfolder << endl;
        int status;
        status = mkdir(checknewoutputPrefixfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (status == -1) cout << "FAILED TO MAKE NEW FOLDER" << endl;
        else cout << "MAKE DIRECTORY" << endl;
    }
    else{
        cout << "LATTICE : " << checknewoutputPrefixfolder << " IS NOT A DIRECTORY" << endl;
        int status;
        status = mkdir(checknewoutputPrefixfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (status == -1) cout << "FAILED TO MAKE NEW FOLDER" << endl;
        else cout << "MAKE DIRECTORY" << endl;
    }
    OutputFile crdFile(newoutputPrefixfolder+"/"+newoutputPrefixfile+"_crds.dat");
    OutputFile rcrdFile(newoutputPrefixfolder+"/"+newoutputPrefixfile+"_rcrds.dat");
    //Set up output files

    crdFile.initVariables(6,3,60,12);
    rcrdFile.initVariables(6,3,60,12);

//    cout << latticeCnd << endl;
//    cout << nodeCnd << endl;

    //Write coordination
    crdFile.write(latticeCnd,nodeCnd);

    //Write periodicity
    crdFile.writeRowVector(periodicity);
    rcrdFile.writeRowVector(periodicity);

    //Write node coordinates
    for(int i=0; i<nodes.n; ++i) {
        crdFile.writeRowVector(nodes[i].crd);
//        cout << "NODE  " << i << endl;
//        cout << "    " << nodes[i].crd[0] << "   " << nodes[i].crd[1] << endl;
    }

    //Write ring coordinates
    for(int i=0; i<rings.n; ++i) rcrdFile.writeRowVector(rings[i].crd);
}


void Lattice::writeNetwork(string prefix, int index) {
    //Write node connections and rings
    int dummyfour = int(nodes.n * fractionH);
    int dummytwo = int(nodes.n * fractionL);
    int dummythree = nodes.n - dummyfour - dummytwo;

    if (dummytwo % 2 != 0) {
        dummytwo += 1;
        dummythree -= 1;
    }
    if (dummyfour% 2 != 0) {
        dummyfour+= 1;
        dummythree-= 1;
    }
    if ((nodes.n - dummyfour- dummythree- dummytwo) == -1) {
        dummythree-= 1;
    }
    else if ((nodes.n - dummyfour- dummythree- dummytwo) == +1) {
        dummythree += 1;
    }

    int dummylinear;
    if (useLlinear == 1) dummylinear = int(dummytwo*fractionLlinear);
    else dummylinear=dummytwo*0.5;


    string newoutputPrefixfolder;
    string newoutputPrefixfile;
    if (usearray == 1){
        newoutputPrefixfolder = string("./output_"+outputfolder);
        newoutputPrefixfile = string(secondaryinputFile+"_out");
    }
    else if(outStyle==1) {
        newoutputPrefixfolder = "./output/"+outputPrefix +"_4_" + to_string(dummyfour)+"_3_"+to_string(dummythree)+"_2_"+to_string(dummytwo)+"_linear_"+to_string(dummylinear);
        newoutputPrefixfile = "4_" + to_string(dummyfour)+"_3_"+to_string(dummythree)+"_2_"+to_string(dummytwo) +"_linear_"+to_string(dummylinear);
    }
    else {
        newoutputPrefixfolder = "./output";
        newoutputPrefixfile = outputPrefix;
    }

    const char *checknewoutputPrefixfolder = newoutputPrefixfolder.c_str();

    struct stat info;
    if (stat(checknewoutputPrefixfolder, &info) !=0){
        cout << "LATTICE : CANNOT ACCESS DIRECTORY " << checknewoutputPrefixfolder << endl;
        cout << checknewoutputPrefixfolder << endl;
        cout << newoutputPrefixfolder << endl;
        cout << newoutputPrefixfile << endl;
        int status;
        status = mkdir(checknewoutputPrefixfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (status == -1) cout << "FAILED TO MAKE NEW FOLDER" << endl;
        else cout << "MAKE DIRECTORY" << endl;
    }
    else {
        cout << "LATTICE : " << checknewoutputPrefixfolder << " IS NOT A DIRECTORY" << endl;
        int status;
        status = mkdir(checknewoutputPrefixfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (status == -1) cout << "FAILED TO MAKE NEW FOLDER" << endl;
        else cout << "MAKE DIRECTORY" << endl;
    }
   //Set up output file

    OutputFile sampleFile(newoutputPrefixfolder+"/"+newoutputPrefixfile+"_sample_"+to_string(index)+".dat");

    //Write node connections 
    sampleFile.write(nodes.n);
    for(int i=0; i<nodes.n; ++i) {
        sampleFile.writeRowVector(nodes[i].cnxs);
    }
    //Write ring ids
    sampleFile.write(Ring::totalActive);
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) sampleFile.write(i);
    }
    //Write ring connections
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) sampleFile.writeRowVector(rings[i].cnxs);
    }
    //Write rings
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) sampleFile.writeRowVector(rings[i].nodes);
    }
    //Write ring centres of mass
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) sampleFile.writeRowVector(rings[i].crd);
    }
    //Write chains
    sampleFile.write(chains.n);
    for(int i=0; i<chains.n; ++i){
        sampleFile.writeRowVector(chains[i].nodes);
    }
    //Write ring environments
    if(calcEnvs) for(int i=0; i<nodes.n; ++i) sampleFile.write(rings[i].environment);



}

void Lattice::Monolayer(string prefix){
    cout << " ##### Starting Monolayer" << endl;
    float si_si_distance = 1.609* sqrt((32.0/9.0));
    cout << "sisi    " << si_si_distance << endl;
//    float si_o_distance = si_si_distance/2.0;
    float o_o_distance =  1.609*sqrt((8.0/3.0));
    float h = sin((19.5/180)*M_PI)*1.609;




    int si_0, si_1;

    VecF<double> diff;
    int oxygen_number = 0;

    int dummyfour = int(nodes.n * fractionH);
    int dummytwo = int(nodes.n * fractionL);
    int dummythree = nodes.n - dummyfour - dummytwo;

    if (dummytwo % 2 != 0) {
        dummytwo += 1;
        dummythree -= 1;
    }
    if (dummyfour% 2 != 0) {
        dummyfour+= 1;
        dummythree-= 1;
    }
    if ((nodes.n - dummyfour- dummythree- dummytwo) == -1) {
        dummythree-= 1;
    }
    else if ((nodes.n - dummyfour- dummythree- dummytwo) == +1) {
        dummythree += 1;
    }

    int dummylinear;
    if (useLlinear == 1) dummylinear = int(dummytwo*fractionLlinear);
    else dummylinear=dummytwo*0.5;

    string newoutputPrefixfolder;
    string newoutputPrefixfile;
    if (usearray == 1){
        newoutputPrefixfolder = string("./output_"+outputfolder);
        newoutputPrefixfile = string(secondaryinputFile+"_out");
    }
    else if(outStyle==1) {
        newoutputPrefixfolder = "./output/"+outputPrefix +"_4_" + to_string(dummyfour)+"_3_"+to_string(dummythree)+"_2_"+to_string(dummytwo)+"_linear_"+to_string(dummylinear);
        newoutputPrefixfile = "4_" + to_string(dummyfour)+"_3_"+to_string(dummythree)+"_2_"+to_string(dummytwo) +"_linear_"+to_string(dummylinear);
    }
    else {
        newoutputPrefixfolder = "./output";
        newoutputPrefixfile = outputPrefix;
    }

    double dim_x, dim_y;
    dim_x = (latticeDim*1.5)*si_si_distance;
    dim_y = (latticeDim*0.5)*sqrt(3.0)*si_si_distance;
    OutputFile dimFile(newoutputPrefixfolder+"/dimensions.dat");
    dimFile.initVariables(8,3,60,15);
    dimFile.write(0.5*latticeDim*latticeDim*1.5*sqrt(3.0)*si_si_distance*si_si_distance);
    dimFile.write(dim_x);
    dimFile.write(dim_y);
    dimFile.write(0.5*latticeDim*latticeDim*1.5*sqrt(3.0)*si_si_distance*si_si_distance / (dim_x*dim_y));
    OutputFile crysFile(newoutputPrefixfolder+"/crys.crds");
    OutputFile harmpairsFile(newoutputPrefixfolder+"/harmpairs.dat");

    VecF<double> Si(3);
    VecF<double> O(3);
    VecF<int> Pair(2);

    int Si_O_harmpairs [2*nodes.n] [5];
    for (int i=0;i<2*nodes.n;++i){
        for (int j=0;j<5;++j){
            Si_O_harmpairs[i][j]=0;
        }
    }

    crysFile.initVariables(8,3,60,15);
    harmpairsFile.initVariables(6,10,60,8);

    harmpairsFile.write(10*nodes.n*2);

    int atom_count = 1;
    cout << "##### Si crys" << endl;
    for (int i=0;i<nodes.n;++i){
        Si_O_harmpairs[2*i][0] = 2*i+1;
        Si_O_harmpairs[2*i+1][0] = 2*i+2;

        Si[0] = nodes[i].crd[0]*si_si_distance;
        Si[1] = nodes[i].crd[1]*si_si_distance;
        Si[2] = 5.0;
        crysFile.writeRowVector(Si);
        atom_count += 1;
        Si[2] = 5.0+2.0*1.609;
        crysFile.writeRowVector(Si);
        atom_count += 1;
    }
    cout << "##### O crys" << endl;
    for (int i=0;i<nodes.n;++i){
        if (nodes[i].cnxs.n==3){
            O[0] = nodes[i].crd[0]*si_si_distance; // needs pbc work xoxo

            if (O[0]<0) O[0] += si_si_distance*latticeDim;
            else if (O[0]>si_si_distance*latticeDim) O[0] -= si_si_distance*latticeDim;

            O[1] = nodes[i].crd[1]*si_si_distance;

            if (O[1]<0) O[1] += si_si_distance*latticeDim;
            else if (O[1]>si_si_distance*latticeDim) O[1] -= si_si_distance*latticeDim;

            O[2] = 5+1.609;
            crysFile.writeRowVector(O);
            Pair[0] = atom_count;
            Pair[1] = 2*i+1;
            harmpairsFile.writeRowVector(Pair);
            Si_O_harmpairs[2*i][1] = atom_count;


            Pair[1] = 2*i+2;
            harmpairsFile.writeRowVector(Pair);
            Si_O_harmpairs[2*i+1][1] = atom_count;
            atom_count += 1;
        }
    }

    oxygen_number = 0;
    for(int i=0; i<nodes.n; ++i) {
        for (int j = 0; j < nodes[i].cnxs.n; ++j) {
            if (nodes[i].cnxs[j] > i) {
                diff = nodes[nodes[i].cnxs[j]].crd - nodes[i].crd;
                float mic = latticeDim / 2;
                if (diff[0] > mic) diff[0] -= latticeDim;
                else if (diff[0] < -mic) diff[0] += latticeDim;
                if (diff[1] > mic) diff[1] -= latticeDim;
                else if (diff[1] < -mic) diff[1] += latticeDim;

                O[0] = (nodes[i].crd[0] + diff[0] / 2.) * si_si_distance;

                if (O[0]<0) O[0] += si_si_distance*latticeDim;
                else if (O[0]>si_si_distance*latticeDim) O[0] -= si_si_distance*latticeDim;

                O[1] = (nodes[i].crd[1] + diff[1] / 2.) * si_si_distance;

                if (O[1]<0) O[1] += si_si_distance*latticeDim;
                else if (O[1]>si_si_distance*latticeDim) O[1] -= si_si_distance*latticeDim;

                O[2] = 5 + 2 * 1.609 + h;

                crysFile.writeRowVector(O);
                Pair[0] = atom_count;
                Pair[1] = 2*i+1;
                harmpairsFile.writeRowVector(Pair);
                Pair[1] = 2*nodes[i].cnxs[j]+1;
                harmpairsFile.writeRowVector(Pair);
                int k=1;
                while (k<5){
                    if (Si_O_harmpairs[2*i][k]==0){
                        Si_O_harmpairs[2*i][k] = atom_count;
                        k +=5;
                    }
                    else ++k;
                }
                k=1;
                while (k<5){
                    if (Si_O_harmpairs[2*nodes[i].cnxs[j]][k]==0){
                        Si_O_harmpairs[2*nodes[i].cnxs[j]][k] = atom_count;
                        k +=5;
                    }
                    else ++k;
                }
                atom_count +=1;

                O[2] = 5 - h;
                crysFile.writeRowVector(O);
                Pair[0] = atom_count;
                Pair[1] = 2*i+2;
                harmpairsFile.writeRowVector(Pair);
                Pair[1] = 2*nodes[i].cnxs[j]+2;
                harmpairsFile.writeRowVector(Pair);
                k=1;
                while (k<5){
                    if (Si_O_harmpairs[2*i+1][k]==0){
                        Si_O_harmpairs[2*i+1][k] = atom_count;
                        k +=5;
                    }
                    else ++k;
                }
                k=1;
                while (k<5){
                    if (Si_O_harmpairs[2*nodes[i].cnxs[j]+1][k]==0){
                        Si_O_harmpairs[2*nodes[i].cnxs[j]+1][k] = atom_count;
                        k +=5;
                    }
                    else ++k;
                }
                atom_count +=1;
            }
        }
    }
    cout << "##### O-O pairs" << endl;
    for (int i=0;i<2*nodes.n;++i){
        Pair[0] = Si_O_harmpairs[i][1];
        Pair[1] = Si_O_harmpairs[i][2];
        harmpairsFile.writeRowVector(Pair);
        Pair[0] = Si_O_harmpairs[i][1];
        Pair[1] = Si_O_harmpairs[i][3];
        harmpairsFile.writeRowVector(Pair);
        Pair[0] = Si_O_harmpairs[i][1];
        Pair[1] = Si_O_harmpairs[i][4];
        harmpairsFile.writeRowVector(Pair);
        Pair[0] = Si_O_harmpairs[i][2];
        Pair[1] = Si_O_harmpairs[i][3];
        harmpairsFile.writeRowVector(Pair);
        Pair[0] = Si_O_harmpairs[i][2];
        Pair[1] = Si_O_harmpairs[i][4];
        harmpairsFile.writeRowVector(Pair);
        Pair[0] = Si_O_harmpairs[i][3];
        Pair[1] = Si_O_harmpairs[i][4];
        harmpairsFile.writeRowVector(Pair);
    }
    cout << "##### Tetrahedron" << endl;
    OutputFile tetrehedraFile(newoutputPrefixfolder+"/tetrahedra.dat");
    tetrehedraFile.initVariables(6,3,60,12);
    VecF<int> Tetrahedron(5);
    for (int i=0;i<2*nodes.n;++i){
        for (int j=0;j<5;j++){
            Tetrahedron[j] = Si_O_harmpairs[i][j];
        }
        tetrehedraFile.writeRowVector(Tetrahedron);
    }
    cout << "##### Rings" << endl;
    OutputFile siringsFile(newoutputPrefixfolder+"/rings_si.dat");
    OutputFile oringsFile(newoutputPrefixfolder+"/rings_o.dat");
    OutputFile sioringsFile(newoutputPrefixfolder+"/rings_si_o.dat");

    siringsFile.initVariables(6,3,60,12);
    oringsFile.initVariables(6,3,60,12);
    sioringsFile.initVariables(6,3,60,12);
    VecR<int> Si_Ring (0,10);
    VecR<int> O_Ring (0,10);
    VecR<int> Si_O_Ring (0,20);

    //Write rings
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) {
            Si_Ring.setSize(rings[i].nodes.n);
            O_Ring.setSize(rings[i].nodes.n);
            Si_O_Ring.setSize(2*rings[i].nodes.n);
            for (int j=0;j<rings[i].nodes.n;++j){
                Si_Ring[j] = 2*rings[i].nodes[j];

                if (j==0) si_0 = rings[i].nodes[rings[i].nodes.n-1];
                else si_0 = rings[i].nodes[j-1];
                si_1 = rings[i].nodes[j];

                Si_O_Ring[2*j] = 2*si_1;
                for (int k0=0;k0<5;++k0){
                    for (int k1=0;k1<5;++k1){
                        if (Si_O_harmpairs[2*si_0][k0]==Si_O_harmpairs[2*si_1][k1]){
                            O_Ring[j] = Si_O_harmpairs[2*si_0][k0];
                            Si_O_Ring[2*j+1] = Si_O_harmpairs[2*si_0][k0];
                        }
                    }
                }
            }
            siringsFile.writeRowVector(Si_Ring);
            sioringsFile.writeRowVector(Si_O_Ring);
            oringsFile.writeRowVector(O_Ring);


            for (int j=0;j<rings[i].nodes.n;++j){
                Si_Ring[j] = 2*rings[i].nodes[j]+1;

                if (j==0) si_0 = rings[i].nodes[rings[i].nodes.n-1];
                else si_0 = rings[i].nodes[j-1];
                si_1 = rings[i].nodes[j];

                Si_O_Ring[2*j] = 2*si_1+1;
                for (int k0=0;k0<5;++k0){
                    for (int k1=0;k1<5;++k1){
                        if (Si_O_harmpairs[2*si_0+1][k0]==Si_O_harmpairs[2*si_1+1][k1]){
                            O_Ring[j] = Si_O_harmpairs[2*si_0+1][k0];
                            Si_O_Ring[2*j+1] = Si_O_harmpairs[2*si_0+1][k0];
                        }
                    }
                }
            }
            siringsFile.writeRowVector(Si_Ring);
            sioringsFile.writeRowVector(Si_O_Ring);
            oringsFile.writeRowVector(O_Ring);

        }

    }
    OutputFile ringcnxsFile(newoutputPrefixfolder+"/rings_connections.dat");
    ringcnxsFile.initVariables(6,3,60,12);
    VecR<int> Ring_Cnxs (0,10);
    //Write ring connections
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) {
            Ring_Cnxs.setSize(rings[i].cnxs.n);
            for (int j=0;j<rings[i].cnxs.n;++j){
                Ring_Cnxs[j] = 2*rings[i].cnxs[j];
            }
            ringcnxsFile.writeRowVector(Ring_Cnxs);
            for (int j=0;j<rings[i].cnxs.n;++j){
                Ring_Cnxs[j] = 2*rings[i].cnxs[j]+1;
            }
            ringcnxsFile.writeRowVector(Ring_Cnxs);
        }
    }

    cout << "##### LJ Pairs" << endl;
    OutputFile ljpairsFile(newoutputPrefixfolder+"/lj_pairs.dat");
    ljpairsFile.initVariables(6,3,60,12);
    for (int i=0;i<nodes.n;++i){
        Pair[0] = 2*i;
        Pair[1] = 2*i+1;
        ljpairsFile.writeRowVector(Pair);
    }
    for (int i=0;i<rings.n;++i){
        if (rings[i].active){
            for (int j=0;j<rings[i].nodes.n;j++){
                for (int k=j;k<rings[i].nodes.n;k++){
                    if (rings[i].nodes[j]>rings[i].nodes[k]) {
                        Pair[0] = 2 * rings[i].nodes[j];
                        Pair[1] = 2 * rings[i].nodes[k];
                        ljpairsFile.writeRowVector(Pair);
                        Pair[0] = 2 * rings[i].nodes[j] + 1;
                        Pair[1] = 2 * rings[i].nodes[k] + 1;
                        ljpairsFile.writeRowVector(Pair);
                    }
                }
            }
        }
    }


    cout << "##### OPTIMISE_SILICA_INPUT" << endl;
    OutputFile optimise_silica_File(newoutputPrefixfolder+"/optimise_silica.inpt");
    optimise_silica_File.initVariables(6,3,60,12);
    //
    optimise_silica_File.write("I/O");
    optimise_silica_File.write("./crys.crds              input coordinates");
    optimise_silica_File.write("./harmpairs.dat              input harmonic pairs");
    optimise_silica_File.write("./lj_pairs.dat              input repulsive pairs");
    optimise_silica_File.write("./fixedz.dat              input fixed z atoms");
    optimise_silica_File.write("./test                              output prefix");
    optimise_silica_File.write("-----------------------------------------------------------");
    optimise_silica_File.write("Restart Options");
    optimise_silica_File.write("0               print restart file");
    optimise_silica_File.write("0               restart run?");
    optimise_silica_File.write("-----------------------------------------------------------");
    optimise_silica_File.write("Sample Information");
    optimise_silica_File.write("800                 number Si atoms");
    optimise_silica_File.write("1               stretch x");
    optimise_silica_File.write("1               stretch y");
    optimise_silica_File.write("0               stretch z");
    optimise_silica_File.write("45              central angle");
    optimise_silica_File.write("0              scan angle");
    optimise_silica_File.write(dim_x);
    optimise_silica_File.write(dim_y);
    optimise_silica_File.write("20              unit cell z");
    optimise_silica_File.write("-----------------------------------------------------------");
    optimise_silica_File.write("Geometry Optimisation");
    optimise_silica_File.write("1                   resize(1/0)");
    optimise_silica_File.write("1.300               starting volume (relative to reference)");
    optimise_silica_File.write("0.900               final volume (relative to reference)");
    optimise_silica_File.write("401                   number of volumes to analyse");
    optimise_silica_File.write("20                  samples per volume");
    optimise_silica_File.write(dim_x);
    optimise_silica_File.write(dim_y);
    optimise_silica_File.write("20                  unit cell z reference");
    optimise_silica_File.write("10000000             max steps iterations steepest descent (per area)");
    optimise_silica_File.write("0.5                 Armijo backtracking constant");
    optimise_silica_File.write("1e-9                convergence tolerance");
    optimise_silica_File.write("-----------------------------------------------------------");
    optimise_silica_File.write("Potential Model");
    optimise_silica_File.write("1                   turn on harmonic interactions (1/0)");
    optimise_silica_File.write("1.609               harmonic Si-O r0");
    optimise_silica_File.write("1                   harmonic Si-O k");
    optimise_silica_File.write("1                   harmonic O-O k");
    optimise_silica_File.write("0                   turn on SiSi harmonics (1/0)");
    optimise_silica_File.write("0                   harmonic Si-Si r0");
    optimise_silica_File.write("0                   harmonic Si-Si k");
    optimise_silica_File.write("0                   turn on 24-12 repulsions (1/0)");
    optimise_silica_File.write("2.4                 repulsive r0");
    optimise_silica_File.write("0.25                repulsive k");
    optimise_silica_File.write("0                   turn on z fixing (1/0)");
    optimise_silica_File.write("1                   turn on multicore optimisation");
    optimise_silica_File.write("0                   Parallelise Sample");
    optimise_silica_File.write("1                   Parallelise Area");
    optimise_silica_File.write("1                   number of cores");
    optimise_silica_File.write("0                   turn on cuda");
    optimise_silica_File.write("-----------------------------------------------------------");
    //


    cout << "FINISHED" << endl;

}

int Lattice::pipeenergytot(int latdim) {
    latticeDim = latdim;
    int e = 0;
    for (int i = 0; i < latticeDim * latticeDim; i++) {
        int x = int(i % latticeDim);
        int y = (int((i - x) / latticeDim));
        pipe &p = pipes[i];

        //

        Vector2i v;

        v.x=x;
        v.y=y;

        int ev = 0;
        for(auto d:SecondDIR) {
            if (pipes[v.x+latticeDim*v.y].isConnect(d)) {
                Vector2i neighbourpos = putin(v+d, latticeDim);
                if (!pipes[neighbourpos.x+latticeDim*neighbourpos.y].isConnect(-d)){
                    ev +=1;
                }
            }
        }
        e += ev;


        //

    }
    return e;
}

int Lattice::pipeenergyVector(Vector2i v, int latdim) {
    int N = latdim;
    int ev = 0;
    for(auto d:SecondDIR) {
        if (pipes[v.x+N*v.y].isConnect(d)) {
//        if (cell(v, grid).isConnect(d)) {
            Vector2i neighbourpos = putin(v+d, N);
            if (!pipes[neighbourpos.x+N*neighbourpos.y].isConnect(-d)){
//            if (!cell(neighbourpos, grid).isConnect(-d)) {
                ev +=1;
            }
        }
    }
    return ev;
}

int Lattice::pipeenergynode(int pipeid, int latdim) {
    int id = pipeid;
    Vector2i v;
    int latticedim = latdim;
    v.x=int(id%latticedim);
    v.y=int((id-v.x)/latticedim);

    int ev = 0;
    for(auto d:SecondDIR) {
//        if (cell(v, grid).isConnect(d)) {
        if (pipes[v.x+latticedim*v.y].isConnect(d)) {
            Vector2i neighbourpos = putin(v+d, latticedim);
//            if (!cell(neighbourpos, grid).isConnect(-d)) {
            if (!pipes[neighbourpos.x+latticedim*neighbourpos.y].isConnect(-d)){
                ev +=1;
            }
        }
    }
    return ev;
}