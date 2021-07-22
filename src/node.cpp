#include "node.h"




int Node::autoId = 0;

Node::Node() {}


Node::Node(int maxCnxs, float fracL) {
    //Construct with maximum allowed nodes

    //Set incremental id
    numconnections = maxCnxs;
    id=Node::autoId;
    ++Node::autoId;
    fractionLlinear = fracL;
//    cout << "NODE FRACTION LINEAR    " << fractionLlinear << endl;
    //Initialise containters
    allowedCnxs=VecR<int>(0,maxCnxs);
//    secallowedCnxs=VecR<int>(0,maxCnxs-1);
    cnxs=VecR<int>(0,maxCnxs);
    rings=VecR<int>(0,maxCnxs);
    crd=VecF<double>(2);
    environment=-1;
    cluster=-1;
}


void Node::resetRings() {
    //Clear ring vector

    rings=VecR<int>(0,allowedCnxs.n);
}


void Node::setCrd(VecF<double> c) {
    //Set coordinate

    crd=c;
}


void Node::addCnx(int cId) {
    //Add allowed connection id

    allowedCnxs.addValue(cId);
}


void Node::addRing(int rId) {
    //Add ring id

    rings.addValue(rId);
}

//void Node::randomCnxs(int numCnxs, mt19937 &mtGen) {
//    //Shuffle available connections and pick given number
//    cout << "ALLOWED CNX    " << allowedCnxs << endl;
//
//    cnxs=vShuffle(allowedCnxs,mtGen);
//    cnxs.setSize(numCnxs);
//}


void Node::setOrientation(int orientation, int latdim){
    if (cnxs.n==3){
        latticedim = latdim;

        int dummy1,dummy2,dummy3,dummy4;
        int dummysum;
        int A,B,C,D;



        dummy1 = allowedCnxs[0] - id;
        dummy2 = allowedCnxs[1] - id;
        dummy3 = allowedCnxs[2] - id;
        dummy4 = allowedCnxs[3] - id;
        dummysum = (dummy1+dummy2+dummy3+dummy4)/latticedim;

        VecR<int>dummy(4);
        dummy[0]=dummy1;
        dummy[1]=dummy2;
        dummy[2]=dummy3;
        dummy[3]=dummy4;
        if (dummysum==0){
//        cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 4 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==latticedim){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==-latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
        else if (dummysum > 0) {
            // 1,2,3,6
            if (dummysum == 1) {
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 6 "  << endl;
                //A
                for (int j=0;j<4;j++){
                    if (dummy[j]==latticedim-1){
                        A=j;
                    }
                    else if (dummy[j]==-latticedim){
                        D=j;
                    }
                    else if (dummy[j]==1){
                        B=j;
                    }
                    else if (dummy[j]==latticedim){
                        C=j;
                    }
                    else {
                        cout << "FAILS" << endl;
                    }
                }
            }
            else if (dummysum == latticedim-1){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 1 "  << endl;
                for (int j=0;j<4;j++){
                    if (dummy[j]==-1){
                        A=j;
                    }
                    else if (dummy[j]==latticedim){
                        D=j;
                    }
                    else if (dummy[j]==1-latticedim){
                        B=j;
                    }
                    else if (dummy[j]==latticedim*(latticedim-1)){
                        C=j;
                    }
                    else {
                        cout << "FAILS" << endl;
                        cout << dummy[j] << endl;
                    }
                }
            }
            else if (dummysum == latticedim){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 2 "  << endl;
                for (int j=0;j<4;j++){
                    if (dummy[j]==-1){
                        A=j;
                    }
                    else if (dummy[j]==latticedim){
                        D=j;
                    }
                    else if (dummy[j]==1){
                        B=j;
                    }
                    else if (dummy[j]==latticedim*(latticedim-1)){
                        C=j;
                    }
                    else {
                        cout << "FAILS" << endl;
                    }
                }

            }
            else if (dummysum == latticedim+1){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 3 "  << endl;
                for (int j=0;j<4;j++){
                    if (dummy[j]==latticedim-1){
                        A=j;
                    }
                    else if (dummy[j]==latticedim){
                        D=j;
                    }
                    else if (dummy[j]==1){
                        B=j;
                    }
                    else if (dummy[j]==latticedim*(latticedim-1)){
                        C=j;
                    }
                    else {
                        cout << "FAILS" << endl;
                    }
                }
            }
        }
        else {
            // 5,7,8,9
            if (dummysum == -1) {
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 5 "  << endl;
                for (int j=0;j<4;j++){
                    if (dummy[j]==-1){
                        A=j;
                    }
                    else if (dummy[j]==latticedim){
                        D=j;
                    }
                    else if (dummy[j]==1-latticedim){
                        B=j;
                    }
                    else if (dummy[j]==-latticedim){
                        C=j;
                    }
                    else {
                        cout << "FAILS" << endl;
                    }
                }
            }
            else if (dummysum == 1-latticedim){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 9 "  << endl;
                for (int j=0;j<4;j++){
                    if (dummy[j]==latticedim-1){
                        A=j;
                    }
                    else if (dummy[j]==-latticedim*(latticedim-1)){
                        D=j;
                    }
                    else if (dummy[j]==1){
                        B=j;
                    }
                    else if (dummy[j]==-latticedim){
                        C=j;
                    }
                    else {
                        cout << "FAILS" << endl;
                        cout << dummy[j] << endl;
                    }
                }
            }
            else if (dummysum == -latticedim){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 8 "  << endl;
                for (int j=0;j<4;j++){
                    if (dummy[j]==-1){
                        A=j;
                    }
                    else if (dummy[j]==-latticedim*(latticedim-1)){
                        D=j;
                    }
                    else if (dummy[j]==1){
                        B=j;
                    }
                    else if (dummy[j]==-latticedim){
                        C=j;
                    }
                    else {
                        cout << "FAILS" << endl;
                    }
                }
            }
            else if (dummysum == -latticedim-1){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 7 "  << endl;
                for (int j=0;j<4;j++){
                    if (dummy[j]==-1){
                        A=j;
                    }
                    else if (dummy[j]==-latticedim*(latticedim-1)){
                        D=j;
                    }
                    else if (dummy[j]==1-latticedim){
                        B=j;
                    }
                    else if (dummy[j]==-latticedim){
                        C=j;
                    }
                    else {
                        cout << "FAILS" << endl;
                    }
                }
            }
        }


        //orientation translation
        //DIRS = U R D L = D B C A
        //4 -> R D L -> B C A
        //3 -> D L U -> C A D
        //2 -> L U R -> A D B
        //1 -> U R D -> D B C

        if (orientation==4){
            cnxs[0] = allowedCnxs[B];
            cnxs[1] = allowedCnxs[C];
            cnxs[2] = allowedCnxs[A];
        }
        else if (orientation==3){
            cnxs[0] = allowedCnxs[C];
            cnxs[1] = allowedCnxs[A];
            cnxs[2] = allowedCnxs[D];
        }
        else if (orientation==2){
            cnxs[0] = allowedCnxs[A];
            cnxs[1] = allowedCnxs[D];
            cnxs[2] = allowedCnxs[B];
        }
        else if (orientation==1){
            cnxs[0] = allowedCnxs[D];
            cnxs[1] = allowedCnxs[B];
            cnxs[2] = allowedCnxs[C];
        }

    }
    else{
        cout << "WRONG NUMBER OF CONNECTIONS !!" << endl;
    }
}

void Node::clockwiseCnxs(int latdim){
    //include pipe reference for rotation
    latticedim = latdim;

    VecR<int>CurrCnxs(cnxs.n);
    CurrCnxs = cnxs;

    int dummy1,dummy2,dummy3,dummy4;
    int dummysum;
    int A,B,C,D;



    dummy1 = allowedCnxs[0] - id;
    dummy2 = allowedCnxs[1] - id;
    dummy3 = allowedCnxs[2] - id;
    dummy4 = allowedCnxs[3] - id;
    dummysum = (dummy1+dummy2+dummy3+dummy4)/latticedim;

    VecR<int>dummy(4);
    dummy[0]=dummy1;
    dummy[1]=dummy2;
    dummy[2]=dummy3;
    dummy[3]=dummy4;
    if (dummysum==0){
        for (int j=0;j<4;j++){
           if (dummy[j]==-1){
               A=j;
           }
           else if (dummy[j]==latticedim){
               D=j;
           }
           else if (dummy[j]==1){
               B=j;
           }
           else if (dummy[j]==-latticedim){
               C=j;
           }
           else {
               cout << "FAILS" << endl;
           }
        }
    }
    else if (dummysum > 0) {
        // 1,2,3,6
        if (dummysum == 1) {
            //A
            for (int j=0;j<4;j++){
                if (dummy[j]==latticedim-1){
                    A=j;
                }
                else if (dummy[j]==-latticedim){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
        else if (dummysum == latticedim-1){
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==latticedim){
                    D=j;
                }
                else if (dummy[j]==1-latticedim){
                    B=j;
                }
                else if (dummy[j]==latticedim*(latticedim-1)){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                    cout << dummy[j] << endl;
                }
            }
        }
        else if (dummysum == latticedim){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 2 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==latticedim){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==latticedim*(latticedim-1)){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }

        }
        else if (dummysum == latticedim+1){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 3 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==latticedim-1){
                    A=j;
                }
                else if (dummy[j]==latticedim){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==latticedim*(latticedim-1)){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
    }
    else {
        // 5,7,8,9
        if (dummysum == -1) {
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 5 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==latticedim){
                    D=j;
                }
                else if (dummy[j]==1-latticedim){
                    B=j;
                }
                else if (dummy[j]==-latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
        else if (dummysum == 1-latticedim){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 9 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==latticedim-1){
                    A=j;
                }
                else if (dummy[j]==-latticedim*(latticedim-1)){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==-latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                    cout << dummy[j] << endl;
                }
            }
        }
        else if (dummysum == -latticedim){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 8 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==-latticedim*(latticedim-1)){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==-latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
        else if (dummysum == -latticedim-1){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 7 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==-latticedim*(latticedim-1)){
                    D=j;
                }
                else if (dummy[j]==1-latticedim){
                    B=j;
                }
                else if (dummy[j]==-latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
    }

///////////////////////////////////////////////////////////////////////

    int cnxA=10, cnxB=10, cnxC=10, cnxD=10;

    for (int i=0;i<cnxs.n;++i){
        if (cnxs[i]==allowedCnxs[A]){
            cnxA=i;
        }
        else if (cnxs[i]==allowedCnxs[B]){
            cnxB=i;
        }
        else if (cnxs[i]==allowedCnxs[C]){
            cnxC=i;
        }
        else if (cnxs[i]==allowedCnxs[D]){
            cnxD=i;
        }
    }
    cout << "A : " << allowedCnxs[A] << "  B : " << allowedCnxs[B] << "  C : " << allowedCnxs[C] << "  D : " << allowedCnxs[D] << endl;
///////////////////////////////////////////////////////////////////////
    cout << cnxA << "  " << cnxB << "  " << cnxC <<  "  "  << cnxD << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>>>>" << endl;
    cout << "ORIGINAL CNXS" << endl;
    for (int j=0;j<cnxs.n;++j){
        cout << cnxs[j] << "  ";
    }
    cout << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    if (cnxA!=10){
        cnxs[cnxA] = allowedCnxs[C];
    }
    if (cnxB!=10){
        cnxs[cnxB] = allowedCnxs[D];
    }
    if (cnxC!=10){
        cnxs[cnxC] = allowedCnxs[B];
    }
    if (cnxD!=10){
        cnxs[cnxD] = allowedCnxs[A];
    }
///////////////////////////////////////////////////////////////////////
    for (int j=0;j<cnxs.n;++j){
        cout << cnxs[j] << "  ";
    }
    cout << endl;

////    VecR<int>NewCnxs(4);
////    for (int i=0;i<4;i++){
////        if (i==A){
////            NewCnxs[i] = cnxs[D];
////        }
////        else if (i==B){
////            NewCnxs[i] = cnxs[C];
////        }
////        else if (i==C){
////            NewCnxs[i] = cnxs[A];
////        }
////        else if (i==D) {
////            NewCnxs[i] = cnxs[B];
////        }
////    }

//    for (int i=0;i<cnxs.n;i++){
//        cnxs[i] = NewCnxs[i];
//    }
//    cnxs.setSize(numconnections);

}
void Node::randomCnxs(int numCnxs, mt19937 &mtGen, int latdim) {
    //Shuffle available connections and pick given number
    //randomly select Linear with probability related to fraction
    //cout << fractionLlinear << endl;
    latticedim = latdim;
    cnxs=vShuffle(allowedCnxs,mtGen);
    //cout << "ALLOWED CNX    " << cnxs[0] << '    ' << cnxs[1] << '    ' << cnxs[2] << '    '<<cnxs[3]<< endl;

    int dummy1,dummy2,dummy3,dummy4;
    int dummysum;
    int A,B,C,D;


//    cout << "######################" << endl;
//    cout << latticedim << endl;
//    cout << endl;

    dummy1 = cnxs[0] - id;
    dummy2 = cnxs[1] - id;
    dummy3 = cnxs[2] - id;
    dummy4 = cnxs[3] - id;
    dummysum = (dummy1+dummy2+dummy3+dummy4)/latticedim;

    VecR<int>dummy(4);
    dummy[0]=dummy1;
    dummy[1]=dummy2;
    dummy[2]=dummy3;
    dummy[3]=dummy4;

    if (dummysum==0){
//        cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 4 "  << endl;
        for (int j=0;j<4;j++){
           if (dummy[j]==-1){
               A=j;
           }
           else if (dummy[j]==latticedim){
               D=j;
           }
           else if (dummy[j]==1){
               B=j;
           }
           else if (dummy[j]==-latticedim){
               C=j;
           }
           else {
               cout << "FAILS" << endl;
           }
        }
    }
    else if (dummysum > 0) {
        // 1,2,3,6
        if (dummysum == 1) {
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 6 "  << endl;
            //A
            for (int j=0;j<4;j++){
                if (dummy[j]==latticedim-1){
                    A=j;
                }
                else if (dummy[j]==-latticedim){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
        else if (dummysum == latticedim-1){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 1 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==latticedim){
                    D=j;
                }
                else if (dummy[j]==1-latticedim){
                    B=j;
                }
                else if (dummy[j]==latticedim*(latticedim-1)){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                    cout << dummy[j] << endl;
                }
            }
        }
        else if (dummysum == latticedim){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 2 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==latticedim){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==latticedim*(latticedim-1)){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }

        }
        else if (dummysum == latticedim+1){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 3 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==latticedim-1){
                    A=j;
                }
                else if (dummy[j]==latticedim){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==latticedim*(latticedim-1)){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
    }
    else {
        // 5,7,8,9
        if (dummysum == -1) {
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 5 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==latticedim){
                    D=j;
                }
                else if (dummy[j]==1-latticedim){
                    B=j;
                }
                else if (dummy[j]==-latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
        else if (dummysum == 1-latticedim){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 9 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==latticedim-1){
                    A=j;
                }
                else if (dummy[j]==-latticedim*(latticedim-1)){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==-latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                    cout << dummy[j] << endl;
                }
            }
        }
        else if (dummysum == -latticedim){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 8 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==-latticedim*(latticedim-1)){
                    D=j;
                }
                else if (dummy[j]==1){
                    B=j;
                }
                else if (dummy[j]==-latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
        else if (dummysum == -latticedim-1){
//            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ 7 "  << endl;
            for (int j=0;j<4;j++){
                if (dummy[j]==-1){
                    A=j;
                }
                else if (dummy[j]==-latticedim*(latticedim-1)){
                    D=j;
                }
                else if (dummy[j]==1-latticedim){
                    B=j;
                }
                else if (dummy[j]==-latticedim){
                    C=j;
                }
                else {
                    cout << "FAILS" << endl;
                }
            }
        }
    }



    if (numCnxs == 2) {

        uniform_real_distribution<double> distribution(0.0, 1.0);
        float mc = distribution(mtGen);

        VecR<int> todelete;
        todelete = VecR<int>(2);


        if (mc <= fractionLlinear) {
//            cout << "IS MEANT TO BE LINEAR" << endl;

            if ((A==0 or A==1) && (B==0 or B==1)) {
                cnxs.setSize(numCnxs);
//                cout << "ALREADY LINEAR" << endl;
            }
            else if ((C==0 or C==1) && (D==0 or D==1)){
                cnxs.setSize(numCnxs);
//                cout << "ALREADY LINEAR"  << endl;
            }
            else {
                if (A==0 or B==0){
                    if (C>D){
                        cnxs.delValue(cnxs[C]);
                        cnxs.delValue(cnxs[D]);
                    }
                    else {
                        cnxs.delValue(cnxs[D]);
                        cnxs.delValue(cnxs[C]);
                    }
                }
                else if (C==0 or D==0) {
                    if (B>A){
                        cnxs.delValue(cnxs[B]);
                        cnxs.delValue(cnxs[A]);
                    }
                    else {
                        cnxs.delValue(cnxs[A]);
                        cnxs.delValue(cnxs[B]);
                    }
                }

//                cout << cnxs << endl;

            }
        }
        else {
//            cout << "SELECTED TO BE NON LINEAR" << endl;

            if ((A==0 or A==1) && (B==0 or B==1)) {
                cnxs.delValue(cnxs[1]);
                cnxs.setSize(numCnxs);
            }
            else if ((C==0 or C==1) && (D==0 or D==1)){
                cnxs.delValue(cnxs[1]);
                cnxs.setSize(numCnxs);
            }
            else {
                cnxs.setSize(numCnxs);
            }
        }
    }
    else{
        cnxs.setSize(numCnxs);
    }

  //              if (cnxs[0] == id+1 or cnxs[0] == id-1){  //cnxs[0] is horizontal from id
  //  //                cout << "HORIZONTAL" << endl;
  //                  int j = int(numconnections);
  //                  int i = 1;
  //                  while (i < j){
  //  //                    cout << "i : "<< i << " j : " << j << endl;
  //
  //                      if (cnxs[0] != cnxs[i] + 2 or cnxs[0] != cnxs[i] - 2) {
  //                                                                          //nonlinear connections
  //                          cnxs.delValue(cnxs[i]);
  //                          j -=1 ;
  //                      }
  //                      i += 1;
  //                  }
  //              }
  //              else {      //cnxs[0] is vertical from id
  //  //                cout << "VERTICAL" << endl;
  //                  int j = int(numconnections);
  //                  int i = 1;
  //                  while (i < j){
  //  //                    cout << "i : "<< i << " j : " << j << endl;
  //
  //                      if (cnxs[i] == id + 1 or cnxs[i] == id -1) {      //nonlinear connections
  //                          cnxs.delValue(cnxs[i]);
  //                          j -=1 ;
  //                      }
  //                      i+=1;
  //  //                    cout << i << endl;
  //                  }
  //              }
  //  //            cout << "CHOSEN  CNX    " << cnxs[0] << "    " << cnxs[1] << endl;
  //          }
  //          else {
  //              int j = int(numconnections);
  //              int i = 1;
  //              while (i < j){
  //                  if (cnxs[0] == cnxs[i]+2 or cnxs[0] == cnxs[i]-2){      //linear connections
  //                      cnxs.delValue(cnxs[i]);
  //                      j -=1 ;
  //                  i+=1;
  //  //                    cnxs.delValue(cnxs[i]);
  //                  }
  //              }
  //          }
  //          cnxs.setSize(numCnxs);
  //      }
  //      else{
  //          cnxs.setSize(numCnxs);
  //      }

}
void Node::randomCnxs(int numCnxs, VecF<int> pattern, mt19937 &mtGen) {
    //Pick connections in orientation of supplied pattern

    int n=allowedCnxs.n;
    uniform_int_distribution<int> randPos(0,n);
    int patternStart=randPos(mtGen);
    int j=0;
    for(int i=0; i<n; ++i){
        if(pattern[i]==1){

            cnxs[j]=allowedCnxs[(patternStart+i)%n];
            ++j;
        }
    }
    cnxs.setSize(numCnxs);
}


void Node::bendorstraight(int size, int threecoordcount, int fourcoordcount, VecR<int> threecoord, VecR<int> fourcoord){
    if (size == 2) { cout << "NOT YET IMPLEMENTED" << endl; }
    else if (size==4) {
        cout << "Increading coordination to 4" << endl;
        numconnections = 4;
        cout << "956" << endl;
        cnxs.resetMaxSize(4);
        cnxs.setSize(4);
        cout << "958" << endl;
        for (int i=0;i<4;i++){
            int d = allowedCnxs[i];
            if (cnxs[0]!=d && cnxs[1]!=d && cnxs[2]!=d){
                cnxs[3]=d;
            }
        }
        cout << "965" << endl;
        fourcoordcount +=1;
        cout << "967" << endl;
        threecoord -=1;
        cout << "969  ;" << fourcoordcount << " from " << fourcoord.n << " with size limit " << fourcoord.nMax << endl;
        fourcoord.resetMaxSize(fourcoordcount);
        fourcoord.setSize(fourcoordcount);
        cout << "971" << endl;
        fourcoord[-1] = id;
        cout << "973" << endl;
        bool reached = false;
        for (int i=0;i<threecoordcount;++i){
            if (reached==true){
                threecoord[i-1] = threecoord[i];
            }
            if (threecoord[i]==id) {
                cout << " NODE TO REMOVE IS IN (ref) : "<< autoId << endl;
                reached = true;
            }
        }
        threecoord.setSize(threecoordcount);
        cout << "975" << endl;
    }
    else if (size==3){
        cout << "Reducing coordination to 3" << endl;
        numconnections = 3;
        cnxs.setSize(3);
        fourcoordcount -=1;
        threecoord +=1;
        fourcoord.delValue(autoId);
        threecoord.setSize(threecoordcount);
        threecoord[-1] = autoId;

    }
    cout << "made it out of bendorstraight" << endl;
}



//

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Node::setCnxs(VecR<int> ids) {
    //Set connections explicitly

    cnxs=ids;
}


void Node::maximiseCnxs() {
    //Set connections to all allowed

    cnxs=allowedCnxs;
}


void Node::breakCnx(int id) {
    //Break connection

    cnxs.delValue(id);
}


void Node::swapRing(int delId, int addId) {
    //Swap all instances of ring id

    for(int i=0; i<rings.n; ++i){
        if(rings[i]==delId) rings[i]=addId;
    }
    rings=vUnique(rings);
}


VecR<int> Node::getReciprocalCnxs(VecR<Node> &others) {
    //Find connections which are reciprocated by other nodes

    VecR<int> rCnxs(0,cnxs.n);
    for(int i=0; i<cnxs.n; ++i){
        if(vContains(others[cnxs[i]].cnxs,id)) rCnxs.addValue(cnxs[i]);
    }
    return rCnxs;
}

VecR<int> Node::getVacantCnxs() {
    //Get allowed connections which are not formed

    VecR<int> vacantCnxs(0, allowedCnxs.n);
    for(int i=0; i<allowedCnxs.n; ++i){
        if(!vContains(cnxs,allowedCnxs[i])) vacantCnxs.addValue(allowedCnxs[i]);
    }
    return vacantCnxs;
}

VecR<int> Node::getSharedRings(Node &other) {
    //Get rings shared with other node

    VecR<int> shared(0,rings.n);
    for(int i=0; i<rings.n; ++i){
        if(vContains(other.rings,rings[i])) shared.addValue(rings[i]);
    }
    return shared;
}

VecR<int> Node::getSharedCnxs(Ring &ring) {
    //Get shared nodes with ring

    VecR<int> shared(0,cnxs.n);
    for(int i=0; i<cnxs.n; ++i){
        if(vContains(ring.nodes,cnxs[i])) shared.addValue(cnxs[i]);
    }
    return shared;
}





sf::Vector2i Up(0,-1);
sf::Vector2i Down(0,1);
sf::Vector2i Right(1,0);
sf::Vector2i Left(-1,0);
sf::Vector2i DIR[4] = {Up,Right,Down,Left};

int pipe::autoId = 0;

pipe::pipe() {}

pipe::pipe(int latdimX, int latdimY) {
    latticedimX = latdimX;
    latticedimY = latdimY;
    angle=0;
    id=pipe::autoId;
    ++pipe::autoId;
    crds.x = int(id%latticedimX);
    crds.y = int((id-crds.x)/latticedimY);

    //replacement for 3 coords
    if ((crds.x+crds.y)%2==0){
        orientation = 1;
        rotate();
        rotate();
    }
    else{
        orientation = 3;
    }
}

void pipe::rotate() {
    for(int i=0;i<dirs.n;i++){
        if (dirs[i]==Up)  dirs[i]=Right;
        else if (dirs[i]==Right) dirs[i]=Down;
        else if (dirs[i]==Down)  dirs[i]=Left;
        else if (dirs[i]==Left)  dirs[i]=Up;
    }
};

Vector2i putin(Vector2i neighbourpos, int latdimX, int latdimY) {
    int nX = latdimX;
    int nY = latdimY;
    //potential neighbour position, given vy v+d
    if (neighbourpos.x<0) neighbourpos.x+=(nX);
    else if (neighbourpos.x>nX-1) neighbourpos.x -=(nX);
    if (neighbourpos.y<0) neighbourpos.y+=(nY);
    else if (neighbourpos.y>nY-1) neighbourpos.y-=(nY);
    return neighbourpos;
}

bool pipe::isConnect(Vector2i dir){
    for(int i=0;i<dirs.n;i++){
        if (dirs[i]==dir) return true;
    }
    return false;
}

void pipe::bendorstraighten() {
    if (dirs.n == 2) {
        std::string s = "";
        for (auto d:DIR) s += isConnect(d) ? '1' : '0';

        if (s == "1010") {
            cout << "       LINEAR1" << endl;
            for (int i = 0; i < 2; i++) {
                if (dirs[i] == Up) {
                    dirs[i] = Right;
                    orientation = 3;
                }
            }
        } else if (s == "0101") {
            cout << "       LINEAR2" << endl;
            for (int i = 0; i < 2; i++) {
                if (dirs[i] == Right) {
                    dirs[i] = Down;
                    orientation = 0;
                }
            }
        } else if (s == "1100") {
            cout << "       Bent1" << endl;
            for (int i = 0; i < 2; i++) {
                if (dirs[i] == Up) {
                    dirs[i] = Left;
                    orientation = 0;
                }
            }
        } else if (s == "0110") {
            cout << "       Bent2" << endl;
            for (int i = 0; i < 2; i++) {
                if (dirs[i] == Right) {
                    dirs[i] = Up;
                    orientation = 1;
                }
            }
        } else if (s == "0011") {
            cout << "       Bent3" << endl;
            for (int i = 0; i < 2; i++) {
                if (dirs[i] == Down) {
                    dirs[i] = Right;
                    orientation = 0;
                }
            }
        } else if (s == "1001") {
            cout << "       Bent4" << endl;
            for (int i = 0; i < 2; i++) {
                if (dirs[i] == Left) {
                    dirs[i] = Down;
                    orientation = 1;
                }
            }
        }

    }
    else if (dirs.n == 3) {
        dirs.resetMaxSize(4);
        dirs.setSize(4);
        if (dirs[0]!=Up && dirs[1]!=Up && dirs[2]!=Up) {dirs[3]=Up;}
        else if (dirs[0]!=Right && dirs[1]!=Right && dirs[2]!=Right) {dirs[3]=Right;}
        else if (dirs[0]!=Left && dirs[1]!=Left && dirs[2]!=Left) {dirs[3]=Left;}
        else if (dirs[0]!=Down && dirs[1]!=Down && dirs[2]!=Down) {dirs[3]=Down;}

        for (int n = 4; n > 0; n--) //find orientation//
        {
            std::string s = "";
            for (auto d: DIR) s += isConnect(d) ? '1' : '0';
            if (s == "0011" || s == "0111" || s == "0101" || s == "1111") orientation = n;
            rotate();
            cout << s << "    ";
        }
        cout << endl;


    }
    else if (dirs.n == 4){
        dirs.setSize(3);
        cout << endl;

        for (int n = 4; n > 0; n--) //find orientation//
        {
            std::string s = "";
            for (auto d: DIR) s += isConnect(d) ? '1' : '0';
            if (s == "0011" || s == "0111" || s == "0101" || s == "1111") orientation = n;
            rotate();
        }
    }
}

void pipe::setCnxs(VecR<int> Nodecnxs) {
    cnxs = VecR<int>(Nodecnxs.n, Nodecnxs.n);
    for (int z=0;z<Nodecnxs.n;z++){

        cnxs[z] = int (Nodecnxs[z]);
    }
}
void pipe::cnxstodirs(){
    int nX, nY;
    nX = latticedimX;
    nY = latticedimY;
//    cout << "Call cnxstodirs" << endl;
    int connection_sum=0;
//    cout << "ID : " << id << endl;
    dirs = VecR<Vector2i>(cnxs.n, cnxs.n);
    for (int z=0; z<cnxs.n; z++){connection_sum+=(cnxs[z]-id)/nX;}
//    cout << "CONNECTION SUM : " << connection_sum << " from " << cnxs.n << " cnxs" << endl;
    for (int z=0; z<cnxs.n; z++){
        int corrected_cnx = cnxs[z]-id;
//        cout << z << "    -> " << corrected_cnx << endl;
        //A
        if (corrected_cnx==-1 || corrected_cnx==nX-1){dirs[z]=Left;}
            //B
        else if (corrected_cnx==1 || corrected_cnx==1-nX){dirs[z]=Right;}
            //C&D
        else if (connection_sum==0) {
//            cout << "sum =0" << endl;
            //C
            if (corrected_cnx == -nX) { dirs[z] = Up; }
                //D
            else if (corrected_cnx == nX) { dirs[z] = Down; }
        }
        else if (connection_sum > 0){
//            cout << "sum >0" << endl;
            if (connection_sum==1) {
                //C
                if (corrected_cnx == -nX) { dirs[z] = Up; }
                    //D
                else if (corrected_cnx == nX) { dirs[z] = Down; }
            }
            else if (connection_sum==nX-1 || connection_sum==nX || connection_sum==nX+1) {
//                cout << "catch 1 " << corrected_cnx/N << endl;
                //C
                if (corrected_cnx == nX*(nX-1)) { dirs[z] = Up; }
                    //D
                else if (corrected_cnx == nX) { dirs[z] = Down; }

//                else {cout << "catch Failed" << endl;}
            }
        }
        else if (connection_sum < 0){
//            cout << "sum <0" << endl;
            if (connection_sum==-1) {
                //C
                if (corrected_cnx == -nX) { dirs[z] = Up; }
                    //D
                else if (corrected_cnx == nX) { dirs[z] = Down; }
            }
            else if (connection_sum==1-nX || connection_sum==-nX || connection_sum==-nX-1) {
//                cout << "catch 1 " << corrected_cnx/N << endl;
                //C
                if (corrected_cnx == -nX*(nX-1)) { dirs[z] = Up; }
                    //D
                else if (corrected_cnx == -nX) { dirs[z] = Down; }
//                else {cout << "catch Failed" << endl;}
            }
        }
//        cout << "( " << dirs[z].x << "," << dirs[z].y << " ) ";
    }
//    cout << "End cnxstodir" << endl;
}

void pipe::setOrientation() {
    for (int n = 4; n > 0; n--) //find orientation//
    {
        std::string s = "";
        for (auto d: DIR) s += isConnect(d) ? '1' : '0';
        if (s == "0011" || s == "0111" || s == "0101" || s == "1111") orientation = n;
        rotate();
    }
}

void pipe::setAngle(){
    angle +=5;
    if (angle > orientation * 90) angle = orientation * 90;
}

int pipe::energynode(VecR<pipe> pipes, int latdimX, int latdimY) {
    Vector2i v;
    latticedimX = latdimX;
    latticedimY = latdimY;
    v.x=int(id%latticedimX);
    v.y=int((id-v.x)/latticedimY);

    int ev = 0;
    for(auto d:DIR) {
//        if (cell(v, grid).isConnect(d)) {
        if (pipes[v.x+latticedimX*v.y].isConnect(d)) {
            Vector2i neighbourpos = putin(v+d, latticedimX, latticedimY);
//            if (!cell(neighbourpos, grid).isConnect(-d)) {
            if (!pipes[neighbourpos.x+latticedimX*neighbourpos.y].isConnect(-d)){
                ev +=1;
            }
        }
    }
    return ev;
}

string pipe::orientationstring(){
    std::string s = "";
    for (auto d: DIR) s += isConnect(d) ? '1' : '0';
//    cout << s << "      Up Right Down Left " << endl;
    return s;
}

int energytot(VecR<pipe> pipes, int latdimX, int latdimY){
    int N = latdim;
    int e = 0;
    for (int i = 0; i < N * N; i++) {
        int x = int(i % N);
        int y = (int((i - x) / N));
        pipe &p = pipes[i];

        //

        Vector2i v;

        v.x=x;
        v.y=y;

        int ev = 0;
        for(auto d:DIR) {
            if (pipes[v.x+N*v.y].isConnect(d)) {
                Vector2i neighbourpos = putin(v+d, N);
                if (!pipes[neighbourpos.x+N*neighbourpos.y].isConnect(-d)){
                    ev +=1;
                }
            }
        }
        e += ev;


        //

    }
    return e;
}

int energyVector(Vector2i v, VecR<pipe> pipes, int latdimX, int latdimY){
    int N = latdim;
    int ev = 0;
    for(auto d:DIR) {
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

vector<int> neighbours(int node, vector<int> display, int latdimX, int latdimY){
    int N = latdim;
    int nodex = int ((node) % N);
    int nodey = (int((node - nodex) / N));
    vector <int> next;
    if (nodex!=N-1) next.push_back(node+1);
    else next.push_back(node-(N-1));

    if (nodex!=0) next.push_back(node-1);
    else next.push_back(node+(N-1));

    if (nodey!=N-1) next.push_back(node+N);
    else next.push_back(nodex);

    if (nodey!=0) next.push_back(node-N);
    else next.push_back((N-1)*N+nodex);
    for (int j=0;j<next.size();j++){
        bool check = true;
        for (int i=0;i<display.size();i++){
            if (next[j]==display[i]) check = false;
        }
        if (check==true){
            display.push_back(next[j]);
        }
    }
    return display;
};

bool isOut(Vector2i v, int latdim) {
    return !IntRect(0, 0, latdim, latdim).contains(v);
}





