/** trailCoreExtension.h */
#ifndef TRAIL_CORES_EXTENSION_H 
#define TRAIL_CORES_EXTENSION_H

#include <map>
#include "trailCore.h"
#include "backwardExtension.h"
#include "forwardExtension.h"

/*
This code uses functions and ideas from KeccakTools
(https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

The functions of this file are used to see if it is possible to extend a 
3-round trail core into a 6-round trail core with weight below some bound. 
*/

class LowWeightExclusion {
    protected:  
        /** The explicitly excluded weights per number of rounds. */
        map<unsigned int, Weight> excludedWeight;
        /* minWeight[nrRounds-1] contains the minimum weight for nrRounds rounds. */
        vector<Weight> minWeight;
    public:
        LowWeightExclusion();
        /** This method tells to exclude trails with weight below the given value
          * for the given number of rounds.
          * @param  nrRounds    The number of rounds.
          * @param  weight  The weight below which trails are to be excluded.
          */
        void excludeBelowWeight(unsigned int nrRounds, const Weight& weight);
        /** For a given number of rounds, this function returns the minimum weight
         * to consider.
         * @param  nrRounds    The number of rounds.
         * @return The minimum weight to consider.
         */
        Weight getMinWeight(unsigned int nrRounds);
        friend ostream& operator<<(ostream& out, const LowWeightExclusion& lwe);
    protected:
        void computeExcludedLowWeight(unsigned int upToNrRounds);
};

/** It looks for all trail cores with @a nrRounds rounds
  * up to total weight @a maxTotalWeight that have the given trail as 
  * suffix or prefix.
  * @param fout       Where to output the found trail cores. 
  * @param trailCore  The starting trail core  
  * @param backwardExtension True if the trail core has to be extended
  *                          in the backward direction, false otherwise.
  * @param nrRounds          The target number of rounds. 
  * @param maxTotalWeight    The maximum total weight for the extended trail core
  * @param minWeightFound    Variable to set the minimum weight of the @nrRounds-round trail cores reached.  
  * @param verbose           If true, the function will display indications about the extension.
  *                     
  *
  */ 
void extendTrailCore(ostream& fout, 
                     const TrailCore& trailCore, 
                     bool backwardExtension,
                     unsigned int nrRounds, 
                     long double maxTotalWeight, 
                     Weight& minWeightFound, 
                     bool verbose);
        
/** This function is like extendTrailCore, except that it 
  * processes all the trails from @a fileNameIn.
  * @param fout        Where to output the found trail cores. 
  * @param fileNameIn  The file that contains the trail core to extend.   
  * @param backwardExtension True if the trail core has to be extended
  *                          in the backward direction, false otherwise.
  * @param nrRounds          The target number of rounds. 
  * @param maxTotalWeight    The maximum total weight for the extended trail core
  * @param minWeightFound    Variable to set the minimum weight of the @nrRounds-round trail cores reached. 
  * 
  */ 
void extendTrailCores(ostream& fout, 
                      const string fileNameIn,
                      bool backwardExtension,
                      unsigned int nrRounds, 
                      long double maxTotalWeight, 
                      Weight& minWeightFound);

void recurseForwardExtendTrailCore(ostream& fout, 
                                   const TrailCore& trailCore, 
                                   unsigned int nrRounds, 
                                   long double maxTotalWeight, 
                                   Weight& minWeightFound, 
                                   bool verbose, 
                                   LowWeightExclusion& knownBounds);
    
void recurseBackwardExtendTrailCore(ostream& fout, 
                                    const TrailCore& trailCore, 
                                    unsigned int nrRounds, 
                                    long double maxTotalWeight, 
                                    Weight& minWeightFound, 
                                    bool verbose, 
                                    LowWeightExclusion& knownBounds);

#endif 
