/** trailCore.cpp */
#include "trailCoreExtension.h"
#include "forwardExtension.h"
#include "trailCore.h"


LowWeightExclusion::LowWeightExclusion()
{
}

void LowWeightExclusion::excludeBelowWeight(unsigned int nrRounds, const Weight& weight)
{
    excludedWeight[nrRounds] = weight;
    minWeight.clear();
}

Weight LowWeightExclusion::getMinWeight(unsigned int nrRounds)
{
    if (nrRounds == 0)
        return 0;
    if (nrRounds > minWeight.size())
        computeExcludedLowWeight(nrRounds);
    return minWeight[nrRounds-1];
}

void LowWeightExclusion::computeExcludedLowWeight(unsigned int upToNrRounds)
{
    minWeight.clear();
    for(unsigned int nrRounds = 1; nrRounds <= upToNrRounds; nrRounds++) {
        map<unsigned int, Weight>::iterator i = excludedWeight.find(nrRounds);
        if (i != excludedWeight.end()) {
            minWeight.push_back(i->second);
        }
        else {
            Weight max = 0;
            for(unsigned int n1=1; n1<=(nrRounds-1); n1++) {
                unsigned int n2 = nrRounds - n1;
                Weight sum = minWeight[n1-1] + minWeight[n2-1];
                if (max < sum) 
                    max = sum;
            }
            minWeight.push_back(max);
        }
    }
}

ostream& operator<<(ostream& out, const LowWeightExclusion& lwe)
{
    for(unsigned int nrRounds = 1; nrRounds <= lwe.minWeight.size(); nrRounds++) {
        map<unsigned int, Weight>::const_iterator i = lwe.excludedWeight.find(nrRounds);
        bool extrapolated = (i == lwe.excludedWeight.end());
        out.width(2); out.fill(' '); out << dec << nrRounds << " rounds: ";
        out.width(3); out.fill(' '); out << lwe.minWeight[nrRounds-1] << " ";
        if (extrapolated)
            out << "+";
        out << endl;
    }
    return out;
}

void extendTrailCore(ostream& fout, const TrailCore& trailCore, 
                     bool backwardExtension, unsigned int nrRounds, 
                     long double maxTotalWeight, Weight& minWeightFound, 
                     bool verbose)
{
    LowWeightExclusion knownBounds; 
    knownBounds.excludeBelowWeight(0, 0);
    knownBounds.excludeBelowWeight(1, 2);
    knownBounds.excludeBelowWeight(2, 8);
    knownBounds.excludeBelowWeight(3, 24);

    if (backwardExtension)
        recurseBackwardExtendTrailCore(fout, trailCore, nrRounds, 
                                        maxTotalWeight, minWeightFound, 
                                        verbose, knownBounds);
    else 
        recurseForwardExtendTrailCore(fout, trailCore, nrRounds, 
                                      maxTotalWeight, minWeightFound, 
                                      verbose, knownBounds);
}

void extendTrailCores(ostream& fout, const string fileNameIn,
                      bool backwardExtension, unsigned int nrRounds, 
                      long double maxTotalWeight, 
                      Weight& minWeightFound)
{
   
    TrailFileIterator trailCores(fileNameIn);
    for (; ! trailCores.isEnd(); ++trailCores) {
        extendTrailCore(fout, *trailCores, backwardExtension, nrRounds,
                          maxTotalWeight, minWeightFound, false);
    }
}

void recurseForwardExtendTrailCore(ostream& fout, 
                                   const TrailCore& trailCore, 
                                   unsigned int nrRounds, 
                                   long double maxTotalWeight, 
                                   Weight& minWeightFound, 
                                   bool verbose, 
                                   LowWeightExclusion& knownBounds)
{
    if (verbose) {
        cout << "recurseExtendTrail Forward ( " ; 
        cout << dec << nrRounds << " rounds, "; 
        cout << "up to weight " << maxTotalWeight; 
        cout << ", from a " << trailCore.nrRounds << "-round trail core of weight ", 
        cout << trailCore.weight << " with wMinDir =  " << trailCore.wMinDir ; 
        cout << " ) "<< endl;
    }
    
    Weight baseWeight = trailCore.weight - trailCore.wMinDir;
    long double maxWeightExtension = maxTotalWeight - (long double) baseWeight
                                     - (long double)knownBounds.getMinWeight(nrRounds - trailCore.nrRounds - 1);
               
    if (maxWeightExtension < (long double) knownBounds.getMinWeight(2)) {
        if (verbose) {
            cout << " --leaving because maxWeightExtension (= " << maxWeightExtension; 
            cout << " ) < knownBounds.getMinWeight(2) (= " ; 
            cout << knownBounds.getMinWeight(2) << " ). "<<  endl; 
        }
        return; 
    }
    ForwardExtensionPreparation prep(trailCore.differences.back(), maxWeightExtension);  
    ForwardExtensionIterator extensions(prep, trailCore.differences.back()); 
    for (; !extensions.isEnd(); ++extensions) {
        ForwardExtension extension = *extensions;
        TrailCore newTrailCore(trailCore, extension);

        // save even if not the desired length yet
        newTrailCore.save(fout);
        if (newTrailCore.nrRounds == nrRounds) { 
            if (newTrailCore.weight < minWeightFound) {
                minWeightFound = newTrailCore.weight; 
                if (verbose)
                    cout << "! " << dec << nrRounds << "-round trail of weight " << newTrailCore.weight << " found" << endl;
            }
            return; 
        } else {
            recurseForwardExtendTrailCore(fout, newTrailCore, nrRounds, maxTotalWeight, minWeightFound, verbose, knownBounds);
        }
    } 
}

void recurseBackwardExtendTrailCore(ostream& fout, const TrailCore& trailCore, 
                                    unsigned int nrRounds, long double maxTotalWeight, 
                                    Weight& minWeightFound, bool verbose, 
                                    LowWeightExclusion& knownBounds)
{
    if (verbose) {
        cout << "recurseExtendTrail Backward ( " ; 
        cout << dec << nrRounds << " rounds, "; 
        cout << "up to weight " << maxTotalWeight; 
        cout << ", from a " << trailCore.nrRounds << "-round trail core of weight ", 
        cout << trailCore.weight << " with wMinRev =  " << trailCore.wMinRev ; 
        cout << " ) "<< endl;
    }
    
    Weight baseWeight = trailCore.weight - trailCore.wMinRev;

    long double maxWeightExtension = maxTotalWeight - (long double) baseWeight
                                     - (long double)knownBounds.getMinWeight(nrRounds - trailCore.nrRounds - 1);
                
    if (maxWeightExtension < (long double) knownBounds.getMinWeight(2)) {
        if (verbose) {
            cout << " --leaving because maxWeightExtension (= " << maxWeightExtension; 
            cout << " ) < knownBounds.getMinWeight(2) (= " ; 
            cout << knownBounds.getMinWeight(2) << " ). "<<  endl; 
        }
        return; 
    }
    BackwardExtensionPreparation prep(trailCore.differences[0], maxWeightExtension);  
    BackwardExtensionIterator extensions(prep, trailCore.differences[0]); 
    for (; !extensions.isEnd(); ++extensions) {
        BackwardExtension extension = *extensions;
        TrailCore newTrailCore(trailCore, extension);
        // save even if not the desired length yet
        newTrailCore.save(fout);
        if (newTrailCore.nrRounds == nrRounds) { 
            if (newTrailCore.weight < minWeightFound) {
                minWeightFound = newTrailCore.weight; 
                if (verbose)
                    cout << "! " << dec << nrRounds << "-round trail of weight " << newTrailCore.weight << " found" << endl;
            }
            return;
        } else {
            recurseBackwardExtendTrailCore(fout, newTrailCore, nrRounds, maxTotalWeight, minWeightFound, verbose, knownBounds);
        }
    } 
}
