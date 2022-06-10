/* 3RoundsTrailCores.cpp */ 
#include "3RoundsTrailCores.h"
#include "KK_trailCores.h"
#include "backwardInKernelExtension.h"
#include "forwardExtension.h"
#include "forwardInKernelExtension.h"
#include "mixedStateIterator.h"
#include "trailCore.h"

enum Parity{N, K};


void checkTrailAndParity(string fileNameIn, vector<Parity> parity)
{
    stringstream stream_tri;
    stream_tri << fileNameIn; 
    stream_tri << "-verif";
    string fileNameOut = stream_tri.str();
 
    ofstream fout(fileNameOut.c_str());
   
    // It removes the trail cores that are 
    set<TrailCore> trailsSet;
    ifstream fin(fileNameIn); 
    TrailFileIterator trails(fileNameIn);
    for (; !trails.isEnd(); ++trails) {
        TrailCore trail = *trails; 
        trail.makeCanonical(); 
        assert(2 * parity.size() == trail.differences.size()); 
        for (unsigned int i = 0; i < parity.size(); i++) {
            if (parity[i] == N && trail.differences[2 * i + 1].isInKernel()) {
                cerr << " This difference shoudl not be in the kernel :  " << endl;
                cerr << trail.differences[2 * i + 1] << endl;
                cerr << " problematic trail :  " << endl;
                cout << trail << endl;  
            }        
            if (parity[i] == K && !trail.differences[2 * i + 1].isInKernel()) {
                cerr << " This difference shoudl be in the kernel :  " << endl;
                cerr << trail.differences[2 * i + 1] << endl;
                cerr << " problematic trail :  " << endl;
                cout << trail << endl; 
            }
        }
        trailsSet.insert(trail);
    } 
    set<TrailCore>::const_iterator it;
    for (it = trailsSet.begin(); it != trailsSet.end(); it++) {
        TrailCore trail = *it;
        (*it).save(fout); 
    }
    fout.close(); 
    produceHumanReadableFile(fileNameOut);
}

KN_TrailCores::KN_TrailCores(int T3, int T1)
:T3(T3), T1(T1) 
{
    stringstream stream_K_TrailCores; 
    stream_K_TrailCores << "K-trailCores-T1-";
    stream_K_TrailCores << T1; 
    stream_K_TrailCores << "-alpha-1-beta-0"; 
    file_K_TrailCores = stream_K_TrailCores.str();
 
    stringstream streamForwardExtension; 
    streamForwardExtension << "KN-trailCores-fromForwardExtension-T3-";
    streamForwardExtension << T3; 
    streamForwardExtension << "-T1-"; 
    streamForwardExtension << T1; 
    fileForwardExtension = streamForwardExtension.str();

    stringstream streamInKernelExtension; 
    streamInKernelExtension << "KN-trailCores-fromInKernelExtension-T3-";
    streamInKernelExtension << T3; 
    streamInKernelExtension << "-T1-"; 
    streamInKernelExtension << T1;
    fileInKernelExtension = streamInKernelExtension.str();
}

void KN_TrailCores::K_TrailCores()
{
    int maxCost; 
    time_t start, ending; 
    ofstream fout(file_K_TrailCores.c_str());  
    time(&start); 

    // The cost is always even
    if (T1 % 2 == 0)
        maxCost = T1; 
    else 
        maxCost = max(0, T1 - 1); 

    traverse_K_TrailCoresTree(maxCost, 1, 0, fout); 
    time(&ending);
    double time = difftime(ending, start);
    fout.close();
    produceHumanReadableFileTwoRoundTrailCores(file_K_TrailCores, 1, 0, true, time);
}

void KN_TrailCores::KN_FromInKernelExtension()
{

    time_t start, ending; 
    time(&start);
    ofstream fout(fileInKernelExtension.c_str());

    unsigned int alpha = 1; 
    unsigned int beta  = 1; 
    // The cost in an even integer less than T3 - T1
    int maxCost2Rounds = T3 - T1 - 1; 
    if (maxCost2Rounds % 2 == 1)
        maxCost2Rounds --;  
    TwoRoundTrailCoreCostBoundFunction costFBareState(alpha, beta); 
    ColumnsSet colSet;
    BareStateCache bareStateCache;

    TroikaState stateC; 
    TroikaState stateD; 

    unsigned int cpt = 0; 
    BareStateIterator bareStates(colSet, bareStateCache, costFBareState, maxCost2Rounds, true);
    for (; !bareStates.isEnd(); ++bareStates) {
        BareState bareState = *bareStates;
        if (bareState.valid) {

            MixedTrailCoreCache mixedCache(bareState);
            ActiveTritsSet activeTritsSet(bareState.firstActiveTritsAllowed); 
            TwoRoundTrailCoreCostFunction costFTrits(alpha, beta);
            N_TrailCore_Iterator iteratorTrits(activeTritsSet, mixedCache,
                                               costFTrits, maxCost2Rounds, false);  
            for (; !iteratorTrits.isEnd(); ++iteratorTrits) {
                TwoRoundTrailCore trailCores = *iteratorTrits;

                long double aMaxWeightExtension = T3 - trailCores.wB;
                BackwardInKernelExtensionPreparation prep(aMaxWeightExtension, trailCores.activeA);

                cpt++; 
                if (cpt % 10000 == 0 )
                    cout << cpt << "-th trail to extend " << endl; 
           
                if (prep.possible) {

                    for (; !trailCores.statesB.isEnd(); ++trailCores.statesB) {
                        stateD = *trailCores.statesB; 
                        stateC.setInvL(stateD);
                        BackwardInKernelExtensionIterator extensions(prep, stateC);
                        for (; !extensions.isEnd(); ++extensions) {
                            BackwardInKernelExtension ext = *extensions;
                            TrailCore trail(ext.stateA, ext.stateB, stateC, stateD, ext.wMinRevA, ext.wBC, trailCores.wB); // TODO changer nom wA et wB pour wMinRev, wMinDir  
                            trail.save(fout);       
                        }
                    }
                }
            }
        }
    }
    time(&ending);
    double time = difftime(ending, start);
    fout.close();
    produceHumanReadableFile(fileInKernelExtension, true, time);
}

void KN_TrailCores::KN_FromForwardExtension()
{
    unsigned int cpt = 0; 
    time_t start, ending; 
    time(&start);
    ofstream fout(fileForwardExtension.c_str());

    try {
        ofstream fout(fileForwardExtension.c_str());
        TrailFileIterator trailsIn(file_K_TrailCores); 
        for (; !trailsIn.isEnd(); ++trailsIn) {
            TrailCore trailToExtend = *trailsIn;
            long double maxWeightExtension = T3 - (long double)trailToExtend.wMinRev; 
            ForwardExtensionPreparation prep(trailToExtend.differences.back(), maxWeightExtension);
            ForwardExtensionIterator extensions(prep, trailToExtend.differences.back());
            cpt++; 
            if (cpt % 1000000 == 0 )
                cout << cpt << "-th trail to extend " << endl; 
            for (; !extensions.isEnd(); ++extensions) {
                ForwardExtension extension = *extensions;
                if (!extension.stateD.isInKernel()) {
                    TrailCore extendedTrail(trailToExtend, extension);
                    extendedTrail.save(fout); 
                }
            } 
        }
    } catch (Exception e) {
        cout << e.reason << endl; 
    }
    time(&ending);
    double time = difftime(ending, start);
    fout.close();
    produceHumanReadableFile(fileForwardExtension, true, time);
}

void KN_TrailCores::nrTrailsFound() 
{
    stringstream fileIn, fileOut; 
    
    fileOut << "KN-" << T3; 
    ofstream fout(fileOut.str());

    fileIn << "KN-trailCores-fromForwardExtension-T3-" << T3 << "-T1-" << T1; 
    TrailFileIterator trails(fileIn.str()); 
    for (; ! trails.isEnd(); ++trails)
        (*trails).save(fout);
    fileIn.clear(); fileIn.str(""); fileIn << "KN-trailCores-fromInKernelExtension-T3-" << T3 << "-T1-" << T1; 
    TrailFileIterator trails_(fileIn.str()); 
    for (; ! trails_.isEnd(); ++trails_)
        (*trails_).save(fout);
    fout.close();
    checkTrailAndParity(fileOut.str(), vector<Parity> {K, N});
}

NK_TrailCores::NK_TrailCores(int T3, int T1):
    T3(T3), T1(T1)
{
    stringstream streamInKernelExtension; 
     streamInKernelExtension << "NK-trailCores-fromInKernelExtension-T3-";
     streamInKernelExtension << T3; 
     streamInKernelExtension << "-T1-"; 
     streamInKernelExtension << T1;
     fileInKernelExtension = streamInKernelExtension.str(); 

     stringstream streamBackwardExtension; 
     streamBackwardExtension << "NK-trailCores-fromBackwardExtension-T3-";
     streamBackwardExtension << T3; 
     streamBackwardExtension << "-T1-"; 
     streamBackwardExtension << T1; 
     fileBackwardExtension = streamBackwardExtension.str();
    
     stringstream stream_K_TrailCores; 
     stream_K_TrailCores << "K-trailCores-T1-";
     stream_K_TrailCores << T1; 
     stream_K_TrailCores << "-alpha-0-beta-1"; 
     file_K_TrailCores = stream_K_TrailCores.str();
}

void NK_TrailCores::K_TrailCores()
{
    int maxCost; 
    time_t start, ending; 
    ofstream fout(file_K_TrailCores.c_str());
    time(&start);

    if (T1 % 2 == 0) // The cost is always even
        maxCost = T1; 
    else 
        maxCost = max(0, T1 - 1); 

    traverse_K_TrailCoresTree(T1, 0, 1, fout); 
    time(&ending);
    double time = difftime(ending, start);
    fout.close();
    produceHumanReadableFileTwoRoundTrailCores(file_K_TrailCores, 0, 1, true, time);
}

void NK_TrailCores::NK_FromInKernelExtension()
{

    time_t start, ending; 
    time(&start);
    ofstream fout(fileInKernelExtension.c_str());

    unsigned int alpha = 1; 
    unsigned int beta  = 1; 
    int maxCost2Rounds = T3 - T1 - 1;
    if (maxCost2Rounds % 2 == 1)
        maxCost2Rounds --;

    TwoRoundTrailCoreCostBoundFunction costFBareState(alpha, beta); 
    ColumnsSet colSet;
    BareStateCache bareStateCache;

    TroikaState stateA; 
    TroikaState stateB; 

    unsigned int cpt = 0; 

    BareStateIterator bareStates(colSet, bareStateCache, costFBareState, maxCost2Rounds, true);
    for (; !bareStates.isEnd(); ++bareStates) {
        BareState bareState = *bareStates;
        if (bareState.valid) {

            MixedTrailCoreCache mixedCache(bareState);
            ActiveTritsSet activeTritsSet(bareState.firstActiveTritsAllowed); 
            TwoRoundTrailCoreCostFunction costFTrits(alpha, beta);
            N_TrailCore_Iterator iteratorTrits(activeTritsSet, mixedCache,
                                               costFTrits, maxCost2Rounds, false);  
            for (; !iteratorTrits.isEnd(); ++iteratorTrits) {
                TwoRoundTrailCore trailCores = *iteratorTrits;

                long double maxWeightExtension = T3 - trailCores.wA;
                ForwardInKernelExtensionPreparation prep(maxWeightExtension, trailCores.activeB);
                cpt ++; 
                if (cpt % 10000 == 0)
                    cout << cpt << "-th trail to extend" << endl; 
                if (prep.couldBeExtended()) {

                    for (; !trailCores.statesB.isEnd(); ++trailCores.statesB) {
                        stateB = *trailCores.statesB; 
                        stateA.setInvL(stateB);
                        TrailCore trailToExtend(stateA, stateB, trailCores.wA, trailCores.wB);


                        ForwardInKernelExtensionIterator extensions(prep, stateB); 
                        for (; !extensions.isEnd(); ++extensions) {
                            ForwardExtension ext = *extensions;
                            TrailCore extendedTrailCore(trailToExtend, ext); 
                            extendedTrailCore.save(fout);
                            
                        }
                    }
                }
            }
        }
    }
    time(&ending);
    cout << cpt << " trails extended" << endl; 
    double time = difftime(ending, start);
    fout.close();
    produceHumanReadableFile(fileInKernelExtension, true, time);
}
 
void NK_TrailCores::NK_FromBackwardExtension()
{
    unsigned int cpt = 0; 
    time_t start, ending; 
    time(&start);
    try {
        TrailFileIterator trailsIn(file_K_TrailCores);
        ofstream fout(fileBackwardExtension.c_str()); 

        for (; !trailsIn.isEnd(); ++trailsIn) {
            TrailCore trailToExtend = *trailsIn; 
            long double maxWeightExtension = T3 - (long double)trailToExtend.wMinDir;
            BackwardExtensionPreparation prep(trailToExtend.differences[0], maxWeightExtension);
            BackwardExtensionIterator extensions(prep, trailToExtend.differences[0]);
            cpt++; 
            if (cpt % 1000000 == 0 )
                cout << cpt << "-th trail to extend" << endl; 
            for (; !extensions.isEnd(); ++extensions) {
                BackwardExtension extension = *extensions;
                if (!extension.stateB.isInKernel()) {
                    TrailCore extendedTrail(trailToExtend, extension);
                    extendedTrail.save(fout);
                }
            }
        }
    } catch(Exception e) {
        cout << e.reason << endl; 
    }
    time(&ending);
    double time = difftime(ending, start);
    produceHumanReadableFile(fileBackwardExtension, true, time);
}

void NK_TrailCores::nrTrailsFound() 
{
    stringstream fileIn, fileOut; 
    
    fileOut << "NK-" << T3; 
    ofstream fout(fileOut.str());

    fileIn << "NK-trailCores-fromBackwardExtension-T3-" << T3 << "-T1-" << T1; 
    TrailFileIterator trails(fileIn.str()); 
    for (; ! trails.isEnd(); ++trails)
        (*trails).save(fout);
    fileIn.clear(); fileIn.str(""); fileIn << "NK-trailCores-fromInKernelExtension-T3-" << T3 << "-T1-" << T1; 
    TrailFileIterator trails_(fileIn.str()); 
    for (; ! trails_.isEnd(); ++trails_)
        (*trails_).save(fout);
    fout.close();
    checkTrailAndParity(fileOut.str(), vector<Parity> {N, K});
}

NN_TrailCores::NN_TrailCores(int T3): T3(T3)
{
    stringstream stream_N_ForForwardExtension; 
    stream_N_ForForwardExtension << "N-trailCores-T3-";
    stream_N_ForForwardExtension << T3; 
    stream_N_ForForwardExtension << "-alpha-2-beta-1.txt";  
    file_N_ForForwardExtension_count = stream_N_ForForwardExtension.str();

    stringstream stream_N_ForBackwardExtension; 
    stream_N_ForBackwardExtension << "N-trailCores-T3-";
    stream_N_ForBackwardExtension << T3 ; 
    stream_N_ForBackwardExtension << "-alpha-1-beta-2.txt";
    file_N_ForBackwardExtension_count = stream_N_ForBackwardExtension.str();

    stringstream streamForwardExtension;
    streamForwardExtension << "NN-trailCores-fromForwardExtension-T3-";
    streamForwardExtension << T3;
    fileForwardExtension = streamForwardExtension.str();

    stringstream streamBackwardExtension;
    streamBackwardExtension << "NN-trailCores-fromBackwardExtension-T3-";
    streamBackwardExtension << T3;
    fileBackwardExtension = streamBackwardExtension.str(); 

}


void NN_TrailCores::NN_FromForwardExtension()
{
    int maxCost; 
    ofstream fout; 
    time_t start, ending; 
    time(&start);

    unsigned int cpt = 0; 
    unsigned int alpha = 2; 
    unsigned int beta = 1; 
    maxCost = T3; 
    if (maxCost % 2 == 1) // The cost is always even
        maxCost--;

    // To count the number of 2-round |N|-trail cores (A, B) of cost 2 wMinRev(A) + wMinDir(B) below maxCost
    UINT64 totalCount_N = 0; 
    UINT64 minCost_N;  
    vector<UINT64> countPerCost_N;


    TroikaState stateB; 
    TroikaState stateA; 

    Sbox sbox;
    TwoRoundTrailCoreCostBoundFunction costFBareState(alpha, beta); 
    ColumnsSet colSet;
    BareStateCache bareStateCache;
    BareStateIterator bareStates(colSet, bareStateCache, costFBareState, maxCost, true);
    

    fout.open(fileForwardExtension.c_str());
    for (; !bareStates.isEnd(); ++bareStates) {
        BareState bareState = *bareStates;
        if (bareState.valid) {

            MixedTrailCoreCache mixedCache(bareState);
            ActiveTritsSet activeTritsSet(bareState.firstActiveTritsAllowed); 
            TwoRoundTrailCoreCostFunction costFTrits(alpha, beta);
            N_TrailCore_Iterator iteratorTrits(activeTritsSet, mixedCache,
                                               costFTrits, maxCost, false);  
            for (; !iteratorTrits.isEnd(); ++iteratorTrits) {
                TwoRoundTrailCore trailCores = *iteratorTrits;   
                // don't save the trail core
                cpt++; 
                if (cpt % 10000 == 0 )
                    cout << cpt << "-th trail to extend " << endl;
                
                for (trailCores.statesB.first(); !trailCores.statesB.isEnd(); ++trailCores.statesB) {
                    stateB = *trailCores.statesB; 
                    stateA.setInvL(stateB); 
                    TrailCore trail(stateA, stateB, trailCores.wA, trailCores.wB); 
                 
                
                    unsigned int cost = alpha * (long double)trail.wMinRev 
                                        + beta * (long double)trail.wMinDir; 
                    if (cost >= countPerCost_N.size())
                            countPerCost_N.resize( cost + 1, 0);
                    countPerCost_N[cost]++; 
                    totalCount_N++; 

                    // forward extension
                    long double maxWeightExtension = T3 - (long double)trail.wMinRev; 
                    ForwardExtensionPreparation prep(stateB, maxWeightExtension);
                    ForwardExtensionIterator extensions(prep, stateB);
                    for (; !extensions.isEnd(); ++extensions) {
                        ForwardExtension extension = *extensions;
                        if (!extension.stateD.isInKernel()) {   
                            TrailCore extendedTrail(trail, extension);
                            extendedTrail.save(fout);             
                        }
                    }                     
                }  
            }
        }
    }
    time(&ending);
    double time = difftime(ending, start);
    fout.close();
    produceHumanReadableFile(fileForwardExtension, true, time);

    fout.open(file_N_ForForwardExtension_count);
    minCost_N = 0;
    while ((minCost_N < countPerCost_N.size()) && (countPerCost_N[minCost_N] == 0))
        minCost_N++;

    fout << dec << totalCount_N << " trails of length 2." << endl;
    
    fout << "Minimum cost: " << dec << minCost_N << endl;
    for (unsigned int i = minCost_N; i < countPerCost_N.size(); i++) {
        if (countPerCost_N[i] > 0) {
            fout.width(8); fout.fill(' ');
            fout << dec << countPerCost_N[i] << " trails of cost ";
            fout.width(2); fout.fill(' ');
            fout << i << endl ;
        }
    }
    fout << endl;
    fout.close();

}

void NN_TrailCores::NN_FromBackwardExtension()
{
    int maxCost;
    ofstream fout; 
    time_t start, ending; 
    time(&start);

    unsigned int alpha = 1; 
    unsigned int beta = 2;
    unsigned int cpt = 0; 
    maxCost = T3 - 1;  // The cost is an even integer below T3 - 1
    if (maxCost % 2 == 1)
        maxCost--;

    // To count the number of 2-round |N|-trail cores (A, B) of cost wMinRev(A) + 2 wMinDir(B) below maxCost
    UINT64 totalCount_N = 0; 
    UINT64 minCost_N;  
    vector<UINT64> countPerCost_N;

    TroikaState stateB; 
    TroikaState stateA; 

    Sbox sbox;
    TwoRoundTrailCoreCostBoundFunction costFBareState(alpha, beta); 
    ColumnsSet colSet;
    BareStateCache bareStateCache;
    BareStateIterator bareStates(colSet, bareStateCache, costFBareState, maxCost, true);
    
    fout.open(fileBackwardExtension.c_str());
    for (; !bareStates.isEnd(); ++bareStates) {
        BareState bareState = *bareStates;
        if (bareState.valid) {

            MixedTrailCoreCache mixedCache(bareState);
            ActiveTritsSet activeTritsSet(bareState.firstActiveTritsAllowed); 
            TwoRoundTrailCoreCostFunction costFTrits(alpha, beta);
            N_TrailCore_Iterator iteratorTrits(activeTritsSet, mixedCache,
                                               costFTrits, maxCost, false);  
            for (; !iteratorTrits.isEnd(); ++iteratorTrits) {
                TwoRoundTrailCore trailCores = *iteratorTrits;   
                // don't save the trail
                cpt++; 
                if (cpt % 10000 == 0 )
                    cout << cpt << "-th trail to extend " << endl;
                
                for (trailCores.statesB.first(); !trailCores.statesB.isEnd(); ++trailCores.statesB) {
                    stateB = *trailCores.statesB; 
                    stateA.setInvL(stateB); 
                    TrailCore trail(stateA, stateB, trailCores.wA, trailCores.wB); 
                    // regarder s'il y un probleme et compter les trail 
                    trail.checkTrailCore(sbox);

                    unsigned int cost = alpha * (long double)trail.wMinRev 
                                        + beta * (long double)trail.wMinDir; 
                    if (cost >= countPerCost_N.size())
                            countPerCost_N.resize( cost + 1, 0);
                    countPerCost_N[cost]++; 
                    totalCount_N++; 

                    // etendre vers la droite
                    long double maxWeightExtension = T3 - (long double)trail.wMinDir; 
                    BackwardExtensionPreparation prep(stateA, maxWeightExtension);
                    BackwardExtensionIterator extensions(prep, stateA);
                    for (; !extensions.isEnd(); ++extensions) {
                        BackwardExtension extension = *extensions;
                        if (!extension.stateB.isInKernel()) {   
                            TrailCore extendedTrail(trail, extension);
                            extendedTrail.save(fout);             
                        }
                    }                     
                }  
            }
        }
    }
    time(&ending);
    double time = difftime(ending, start);
    fout.close();
    produceHumanReadableFile(fileBackwardExtension, true, time);

    fout.open(file_N_ForBackwardExtension_count);
    minCost_N = 0;
    while ((minCost_N < countPerCost_N.size()) && (countPerCost_N[minCost_N] == 0))
        minCost_N++;

    fout << dec << totalCount_N << " trails of length 2." << endl;
    
    fout << "Minimum cost: " << dec << minCost_N << endl;
    for (unsigned int i = minCost_N; i < countPerCost_N.size(); i++) {
        if (countPerCost_N[i] > 0) {
            fout.width(8); fout.fill(' ');
            fout << dec << countPerCost_N[i] << " trails of cost ";
            fout.width(2); fout.fill(' ');
            fout << i << endl ;
        }
    }
    fout << endl;
    fout.close();
}
 
void NN_TrailCores::nrTrailsFound() 
{
    stringstream fileIn, fileOut; 
    
    fileOut << "NN-" << T3; 
    ofstream fout(fileOut.str());

    fileIn << "NN-trailCores-fromForwardExtension-T3-" << T3;
    TrailFileIterator trails(fileIn.str()); 
    for (; ! trails.isEnd(); ++trails)
        (*trails).save(fout);
    fileIn.clear(); fileIn.str(""); fileIn << "NN-trailCores-fromBackwardExtension-T3-" << T3; 
    TrailFileIterator trails_(fileIn.str()); 
    for (; ! trails_.isEnd(); ++trails_)
        (*trails_).save(fout);
    fout.close();
    checkTrailAndParity(fileOut.str(), vector<Parity> {N, N});
}
