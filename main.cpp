/* main.cpp */
#include <fstream>
#include <iostream>
#include <ostream>
#include <time.h>
#include <set>
#include "3RoundsTrailCores.h"
#include "KK_trailCores.h"
#include "trailCoreExtension.h"

int main() {


    /* Generate 3-round trail cores up to weight 35
     * The total number of trail is written in the file : 
     * 1) KK-trailCores-T3-35.txt
     * 1) KN-35-verif.txt
     * 2) NK-35-verif.txt 
     * 3) NN-35-verif.txt 
     * 
     */ 
    
     unsigned int T3 = 35;
     unsigned int T1 = 11; 

    // KK TRAIL CORES
    KK_TrailCores KK(T3); 
    KK.generate_KK_trailCores(); 
    
    // KN TRAIL CORES 
    /*
    KN_TrailCores KN(T3, T1); 
    KN.K_TrailCores(); 
    KN.KN_FromForwardExtension(); 
    KN.KN_FromInKernelExtension(); 
    KN.nrTrailsFound(); 
    */
    
    // NK TRAIL CORES 
    /*
    NK_TrailCores NK(T3, T1); 
    NK.K_TrailCores(); 
    NK.NK_FromBackwardExtension(); 
    NK.NK_FromInKernelExtension(); 
    NK.nrTrailsFound();
    */
   
    // NN TRAIL CORES
    /*
    NN_TrailCores NN(T3); 
    NN.NN_FromBackwardExtension(); 
    NN.NN_FromForwardExtension(); 
    NN.nrTrailsFound();
    */ 

    return 0; 
}

