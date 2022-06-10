/* 3RoundsTrailCores.h */
#ifndef THREE_ROUNDS_TRAIL_CORES_H 
#define THREE_ROUNDS_TRAIL_CORES_H 

/*

This code uses functions and ideas from KeccakTools
(https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

The functions of this file are used to generate all 3-round trail cores of 
parity profile |K|N|, |N|K| or |N|N| up to some weight.

*/
#include "types.h"
#include "bareStateIterator.h"
#include "mixedStateIterator.h"
#include "troikaStateIterator.h"
#include "traversal.h"
#include "trailCore.h"
#include "forwardInKernelExtension.h"
#include "forwardExtension.h"
#include "backwardExtension.h"
#include "backwardInKernelExtension.h"

class KN_TrailCores
{
    private:
        int T3; 
        int T1;
        string fileInKernelExtension; 
        string fileForwardExtension;
        string file_K_TrailCores; 
    public:
        KN_TrailCores(int T3, int T1);
        void K_TrailCores(); 
        void KN_FromInKernelExtension();
        void KN_FromForwardExtension();
        void nrTrailsFound();
};

class NK_TrailCores
{
    private:
        int T3; 
        int T1; 
        string fileInKernelExtension; 
        string fileBackwardExtension;
        string file_K_TrailCores; 
    public:
        NK_TrailCores(int T3, int T1);
        void K_TrailCores(); 
        void NK_FromInKernelExtension();
        void NK_FromBackwardExtension();
        void nrTrailsFound();
};

class NN_TrailCores
{
    private:
        int T3;
        string file_N_ForForwardExtension_count; 
        string file_N_ForBackwardExtension_count;
        string fileForwardExtension; 
        string fileBackwardExtension;  
    public:
        NN_TrailCores(int T3);
        void NN_FromForwardExtension(); 
        void NN_FromBackwardExtension(); 
        void nrTrailsFound();
};

#endif
