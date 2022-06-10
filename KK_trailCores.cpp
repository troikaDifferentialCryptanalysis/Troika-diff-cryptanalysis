/* KK_trailCores.cpp */ 
#include "KK_trailCores.h"
#include "backwardInKernelExtension.h"
#include "state.h"
#include "trailCore.h"
#include "traversal.h"
#include "troikaStateIterator.h"
#include <ctime>
 
KK_TrailCores::KK_TrailCores(long double aT3): T3(aT3)
{
    stringstream stream_KK_TrailCores; 
    stream_KK_TrailCores << "KK-trailCores-T3-";
    stream_KK_TrailCores << (long double) T3; 
    file_KK_TrailCores = stream_KK_TrailCores.str();
}
 
ActiveTritAtCAndD::ActiveTritAtCAndD()
{
    // initialize with the first trit to put in the unit-list 
    type = startingTrit; 
    pC.set(0, 0, 0);
    pD = pC.getSRSL(); 
    yOffset = 0; 
    xOffset = 0; 
    firstPosition = pC; 
}; 

bool ActiveTritAtCAndD::incrementOffSets()
{
    assert(type != startingTrit);
    if (type == onSameTryteColumnAtC
        || type == onTheSameTryteColumnAtCAsTheLastStartingTrit) {

        if (xOffset < 2) {
            xOffset++;  
        } else if (yOffset < 1) { 
            xOffset = 0; 
            yOffset++; 
        } else {
            return false; 
        }
    }

    if (type == onSameColumnAtD 
        || type == onTheSameColumnAtDAsTheLastStartingTrit) {

        if (yOffset < 1) {
            yOffset++;  
        } else 
            return false; 
    }
    setCoordinates();
    return true; 
}

void ActiveTritAtCAndD::setCoordinates()
{

    if (type == startingTrit) {
        pC = firstPosition; 
        pD = pC.getSRSL(); 
    }
    else if (type == onSameTryteColumnAtC
        || type == onTheSameTryteColumnAtCAsTheLastStartingTrit) {
        pC = firstPosition;
        pC.xTranslate(xOffset);
        pC.yTranslate(yOffset); 
        pD = pC.getSRSL();
    }
    else {
        pD = firstPosition; 
        pD.xTranslate(xOffset);
        pD.yTranslate(yOffset); 
        pC = pD.getInvSRSL(); 
    }
}

bool ActiveTritAtCAndD::operator < (const ActiveTritAtCAndD& other) const
{
    if (pC.x < other.pC.x) 
        return true; 
    else if (pC.x == other.pC.x) {
        if (pC.y < other.pC.y)
            return true; 
        else if (pC.y == other.pC.y) {
            if (pC.z < other.pC.z)
                return true; 
        }
    }
    return false; 
}

ostream & operator << (ostream &fout, const ActiveTritAtCAndD &trit)
{ 
    fout << trit.pC << " - " << trit.pD;
    fout << "type : "  << trit.type;
    return fout; 
}

ActiveTritAtCAndD ActiveTritsAtCAndDSet::getCandidateForFirstChildUnit(
                            const vector<ActiveTritAtCAndD>& unitList, 
                            const ActiveStatesCAndDCache& cache) const
{
    ActiveTritAtCAndD child;
    const ActiveTritAtCAndD& parent = unitList.back(); 

    child.type = getFirstChildType(parent, cache);
    
    if (child.type == onSameTryteColumnAtC) {
        child.firstPosition = parent.pC; 
        child.firstPosition.x = 3 * (child.firstPosition.x/ 3);
        child.firstPosition.yTranslate(1);
    }

    if (child.type == onTheSameTryteColumnAtCAsTheLastStartingTrit) {
        child.firstPosition = cache.startingTrits.back().pC; 
        child.firstPosition.x = 3 * (child.firstPosition.x/ 3);
        child.firstPosition.yTranslate(1);
    }

    if (child.type == onSameColumnAtD) {
        child.firstPosition = parent.pD; 
        child.firstPosition.yTranslate(1);
    }

    if (child.type == onTheSameColumnAtDAsTheLastStartingTrit) {
        child.firstPosition = cache.startingTrits.back().pD; 
        child.firstPosition.yTranslate(1);
    }

    if (child.type == startingTrit) {
        // The positions are sorted according to the lexicographic order 
        // [x, y, z]
        TritPosition pC = cache.startingTrits.back().pC; 
        do {
            if (pC.z < 26) {
                (pC.z)++; 
            } else if (pC.y < 2) {
                (pC.z) = 0; 
                (pC.y)++; 
            } else if (pC.x < 8) {
                pC.z = 0; 
                pC.y = 0; 
                (pC.x)++; 
            } else {
                throw EndOfSet();
            }
        } while (cache.stateC.isTritActive(pC));  
        child.firstPosition = pC;
    }
    child.setCoordinates();
    return child; 
}

ActiveTritAtCAndD ActiveTritsAtCAndDSet::getFirstChildUnit(
                            const vector<ActiveTritAtCAndD>& unitList, 
                            const ActiveStatesCAndDCache& cache) const
{
    if (unitList.empty()) 
        return ActiveTritAtCAndD(); 
    ActiveTritAtCAndD child; 
    try {
        child = getCandidateForFirstChildUnit(unitList, cache);
        if (isAValidUnit(child, unitList, cache) == false)
            iterateUnit(unitList, child, cache);
    }
    catch (EndOfSet) {
        throw EndOfSet();
    }
    return child; 
}

void ActiveTritsAtCAndDSet::iterateUnit(
                            const vector<ActiveTritAtCAndD>& unitList, 
                            ActiveTritAtCAndD& current,
                            const ActiveStatesCAndDCache& cache) const
{
    if (current.type != startingTrit) {
        do {
            if ( !current.incrementOffSets())
                throw EndOfSet(); 
        } while (! isAValidUnit(current, unitList, cache));
    }

    else {
        if (unitList.empty()) {
            // The first ActiveTritAtCAndD of the tritsList is being iterated :
            // its position at D is restricted to the first slice
            if (current.firstPosition.setNextXY() == false) 
                throw EndOfSet(); 
            current.setCoordinates(); 
        }
        else {
            // Increments the position according to the lexicographic order 
            // [x, y, z]
            TritPosition& pC = current.firstPosition;
            do {
                if (pC.z < 26) {
                    (pC.z)++; 
                } else if (pC.y < 2) {
                    (pC.z) = 0; 
                    (pC.y)++; 
                } else if (pC.x < 8) {
                    pC.z = 0; 
                    pC.y = 0; 
                    (pC.x)++; 
                } else {
                    throw EndOfSet(); 
                }
            } while (cache.stateC.isTritActive(pC)); 
        }
    }
    current.setCoordinates(); 
}

bool ActiveTritsAtCAndDSet::isAValidUnit(const ActiveTritAtCAndD& trit, 
                                const vector<ActiveTritAtCAndD>& unitList, 
                                const ActiveStatesCAndDCache& cache) const
{ 
    if (unitList.empty()) {
        assert (trit.type == startingTrit); 
        if (trit.pC.z != 0)
            return false; 
        return true; 
    }

    if (trit.type == startingTrit) {
        if (cache.stateC.isTritActive(trit.pC))
            return false;
    } 
    if (trit < cache.startingTrits.back())
        return false; 
    return true; 
}

bool ActiveTritsAtCAndDSet::isCanonical(
                         const vector<ActiveTritAtCAndD>& unitList, 
                         const ActiveStatesCAndDCache& cache) const
{
    if (cache.valid == false)
        return true;
    if (cache.newValidPattern == false)
        return false;
    return true; 
}

enum Type ActiveTritsAtCAndDSet::getFirstChildType( 
                          const ActiveTritAtCAndD& parent, 
                          const ActiveStatesCAndDCache& cache) const
{
    if (cache.hasNeighborAtC(parent) == false)
       return onSameTryteColumnAtC; 
    if (cache.hasNeighborAtD(parent) == false)
        return onSameColumnAtD;
    if (cache.hasNeighborAtC(cache.startingTrits.back()) == false)
        return onTheSameTryteColumnAtCAsTheLastStartingTrit; 
    if (cache.hasNeighborAtD(cache.startingTrits.back()) == false)
        return onTheSameColumnAtDAsTheLastStartingTrit; 
    return startingTrit; 
}

ActiveStatesCAndDCache::ActiveStatesCAndDCache()
: lowestNrActiveTrytesA(0), nrActiveTrytesC(0), nrActiveTrytesD(0),
  dummy(true), valid(false), newValidPattern(false)
{
    for (int x = 0; x < COLUMNS; x++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            for (unsigned int z = 0; z < SLICES; z++) {
                nrTritsOnSameColumnAtD[x][y][z] = 0; 
                nrTritsOnSameTryteColumnAtC[x][y][z] = 0; 
            }
        }
    }
}

void ActiveStatesCAndDCache::pushDummy()
{
    dummy = true; 
}

void ActiveStatesCAndDCache::push(const ActiveTritAtCAndD& trit)
{ 
    dummy = false;

    if (false == stateC.isTheTritInAnActiveTryte(trit.pC))
        nrActiveTrytesC++; 
    if (false == stateD.isTheTritInAnActiveTryte(trit.pD))
        nrActiveTrytesD++; 
    stateC.activateTrit(trit.pC); 
    stateD.activateTrit(trit.pD);

    TritPosition pA; 
    TrytePosition tryteC(trit.pC);
    bool activateNecessarilyTrytesAtA = true; 
    for (unsigned int tritIndex = 0; tritIndex < 3; tritIndex++) {
        pA.set(tryteC, tritIndex);
        pA.invSRSL(); 
        if (possibleTrytesA.isTheTritInAnActiveTryte(pA)) {
            activateNecessarilyTrytesAtA = false;
            break; 
        } 
    }
    if (activateNecessarilyTrytesAtA == true) {
        lowestNrActiveTrytesA++; 
        for (unsigned int tritIndex = 0; tritIndex < 3; tritIndex++) {
            pA.set(tryteC, tritIndex);
            pA.invSRSL();
            possibleTrytesA.activateTrit(pA);
        }
    }


    for (unsigned int y = 0; y < ROWS; y++) {
        if (y != trit.pD.y) {
            nrTritsOnSameColumnAtD[trit.pD.x][y][trit.pD.z] += 1;
        }                  
    } 

    unsigned int x = 3 * (trit.pC.x / 3);
    for (unsigned int y = 0; y < ROWS; y++) {
        if (y != trit.pC.y) {
            for (unsigned int xOffset = 0; xOffset < 3; xOffset++)
                nrTritsOnSameTryteColumnAtC[x + xOffset][y][trit.pC.z] += 1; 
        }
    }
    if (trit.type == startingTrit) 
        startingTrits.push_back(trit);

    // Is it a new valid pattern ?
    valid = false;
    newValidPattern = false; 
    if (isAValidPattern(trit)) {
        ActiveState biggestRepresentative; 
        stateC.setTheBiggestRepresentative(biggestRepresentative);
        valid = true; 
        if (patternsC.find(biggestRepresentative) == patternsC.end()) {
            patternsC.insert(biggestRepresentative);
            newValidPattern = true; 
        } 
    } 
}

void ActiveStatesCAndDCache::pop(const ActiveTritAtCAndD& trit)
{
    if (dummy == true) {
        dummy = false;  
        return ; 
    }

    stateC.deactivateTrit(trit.pC); 
    stateD.deactivateTrit(trit.pD);
    if (false == stateC.isTheTritInAnActiveTryte(trit.pC))
        nrActiveTrytesC--;  
    if (false == stateD.isTheTritInAnActiveTryte(trit.pD))
        nrActiveTrytesD--;
  
    if (possibleTrytesA.isTritActive(trit.pC.getInvSRSL()))
    {
        lowestNrActiveTrytesA--; 
        TritPosition pA; 
        TrytePosition tryteC(trit.pC);
        for (unsigned int tritIndex = 0; tritIndex < 3; tritIndex++) {
            pA.set(tryteC, tritIndex);
            pA.invSRSL(); 
            possibleTrytesA.deactivateTrit(pA); 
        }
    }

    for (unsigned int y = 0; y < ROWS; y++) {
        if (y != trit.pD.y)
            nrTritsOnSameColumnAtD[trit.pD.x][y][trit.pD.z] -= 1;            
    } 

    unsigned int x = 3 * (trit.pC.x / 3);
    for (unsigned int y = 0; y < ROWS; y++) {
        if (y != trit.pC.y) {
            for (unsigned int xOffset = 0; xOffset < 3; xOffset++)
                nrTritsOnSameTryteColumnAtC[x + xOffset][y][trit.pC.z] -= 1; 
        }
    }
    if (trit.type == startingTrit) 
        startingTrits.pop_back();

    valid = false;
    newValidPattern = false; 
}

bool ActiveStatesCAndDCache::hasNeighborAtC(const ActiveTritAtCAndD& trit) const
{
    return (bool) nrTritsOnSameTryteColumnAtC[trit.pC.x][trit.pC.y][trit.pC.z];
}

bool ActiveStatesCAndDCache::hasNeighborAtD(const ActiveTritAtCAndD& trit) const
{
    return (bool) nrTritsOnSameColumnAtD[trit.pD.x][trit.pD.y][trit.pD.z];
}

bool ActiveStatesCAndDCache::isAValidPattern(const ActiveTritAtCAndD& last) const
{

    if (!nrTritsOnSameTryteColumnAtC[last.pC.x][last.pC.y][last.pC.z]) 
        return false; 
    if (!nrTritsOnSameColumnAtD[last.pD.x][last.pD.y][last.pD.z]) 
        return false;
    ActiveTritAtCAndD start = startingTrits.back(); 
    if (!nrTritsOnSameColumnAtD[start.pD.x][start.pD.y][start.pD.z]) 
        return false;   
    if (!nrTritsOnSameTryteColumnAtC[start.pC.x][start.pC.y][start.pC.z]) 
        return false;
    return true; 
}

unsigned int KK_TrailCoreCostFunction::getCost(
                             const vector<ActiveTritAtCAndD>& unitList,
                             const ActiveStatesCAndDCache& cache) const
{
    return 2 * (cache.lowestNrActiveTrytesA
               + cache.nrActiveTrytesC
               + cache.nrActiveTrytesD);
}
 
void ActiveStatesCAndD::set(const vector<ActiveTritAtCAndD>& unitList, 
                            const ActiveStatesCAndDCache& cache, 
                            const KK_TrailCoreCostFunction& costF, 
                            unsigned int aMaxCost)
{
   activeD = cache.stateD;
   activeC = cache.stateC; 
   wMinDirD = 2 * cache.nrActiveTrytesD; 
   wMinRevC = 2 * cache.nrActiveTrytesC; 
}

void KK_TrailCores::generate_KK_trailCores()
{
    time_t start, ending;  
    time(&start); 
    ofstream fout(file_KK_TrailCores.c_str());

    ActiveTritsAtCAndDSet KKSet;
    ActiveStatesCAndDCache KKCache;
    KK_TrailCoreCostFunction KKCostF;
    ActiveStatesCAndDIterator it(KKSet, KKCache, KKCostF, T3, true);

    TroikaState stateC; 
    TroikaState stateD; 

    for (; !it.isEnd(); ++it) {
      
        if (it.cache.valid) {
        
            ActiveStatesCAndD current = *it; 
            TroikaStateIterator statesD(current.activeD);
            long double maxWeightExtension = T3 - current.wMinDirD; 
            BackwardInKernelExtensionPreparation prep(maxWeightExtension, current.activeC);
            
            for (; !statesD.isEnd(); ++statesD) {
                   
                   stateD = *statesD; 
                   stateC.setInvSRSL(stateD);
                   BackwardInKernelExtensionIterator extensions(prep, stateC);
                   for (; !extensions.isEnd(); ++extensions) {
                       BackwardInKernelExtension ext = *extensions;
                       TrailCore trail (ext.stateA, ext.stateB, stateC, stateD, ext.wMinRevA, ext.wBC, current.wMinDirD);
                       trail.save(fout); 
                    }
               }
          }
    }
    time(&ending);
    fout.close();
    double time = difftime(ending, start);
    produceHumanReadableFile(file_KK_TrailCores, true, time);    
}
