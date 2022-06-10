/** mixedStateIterator.cpp */
#include "mixedStateIterator.h"
#include "bareStateIterator.h"
#include "state.h"
#include "trailCore.h"
#include "troikaStateIterator.h"

bool ActiveTrits::incrementCoordinates(const ActiveTrits& firstUnit, bool optimization)
{
    if (x < COLUMNS - 1) 
         x++; 
    else if (z < SLICES - 1) {
        z++;
        if (optimization)
            x = firstUnit.x; 
        else 
            x = 0; 
    } else {
        return false;      
    }
    return true;
}

ActiveTrits ActiveTritsSet::getFirstChildUnit(
                            const vector<ActiveTrits>& unitList, 
                            const ActiveTrailCoreCache& cache) const
{
    (void)cache; 
    ActiveTrits child;
    
    if (unitList.empty()) {
        child.x = 0; 
        child.z = 0; 
        child.isYiActive = firstActiveTritsAllowed[child.getXplus9Z()];
        while (child.isYiActive == 0) {
            if (child.incrementCoordinates(child, false) == false)
                throw EndOfSet(); 
            child.isYiActive = firstActiveTritsAllowed[child.getXplus9Z()];
        }
    } 
    else {
        const ActiveTrits& parent = unitList.back();
        child.x = parent.x; 
        child.z = parent.z;
        // See if we can add an active trit at y = 2.
        if (((parent.isYiActive >> 2) & 0x1) == 0) {
            child.isYiActive = 0x4; 
        }
        // See if we can add trits to the next column 
        else {
            do {
                if (child.incrementCoordinates(unitList[0], cache.kernel) == false)
                    throw EndOfSet(); 
                child.isYiActive = firstActiveTritsAllowed[child.getXplus9Z()];
            } while (child.isYiActive == 0); 
        }
    } 
    return child;
}

ostream & operator << (ostream& fout, const ActiveTrits& trits)
{
    fout << (ColumnPosition) trits; 
    if (trits.isYiActive == 0x3)
        fout << "|X|X| |";
    if (trits.isYiActive == 0x5)
        fout << "|X| |X|";
    if (trits.isYiActive == 0x6)
        fout << "| |X|X|";
    if (trits.isYiActive == 0x2)
        fout << "| |X| |";
    if (trits.isYiActive == 0x4)
        fout << "| | |X |";
    return fout; 
}

void ActiveTritsSet::iterateUnit(const vector <ActiveTrits>& unitList,
                                 ActiveTrits& current,
                                 const ActiveTrailCoreCache& cache) const
{
    (void)unitList; 
    (void)cache;
    if (current.isYiActive == 0x3) 
        current.isYiActive = 0x5; 
    else if (current.isYiActive == 0x5) 
        current.isYiActive = 0x6;  
    else if (current.isYiActive == 0x2) 
        current.isYiActive = 0x4; 
    else {
        do {
            if (current.incrementCoordinates(unitList[0], cache.kernel) == false)
                throw EndOfSet(); 
            current.isYiActive = firstActiveTritsAllowed[current.getXplus9Z()];
        } while (current.isYiActive == 0); 
    }
}

bool ActiveTritsSet::isCanonical(const vector<ActiveTrits>& unitList,
                                 ActiveTrailCoreCache& cache) const
{

    if (cache.kernel == false) 
        return true;

    if (unitList[0].z != 0)
        return false;

    unsigned int lastZ = 0; 
    unsigned int z;

    for (unsigned int i = 0; i < unitList.size(); i++) {
        z = unitList[i].z; 
        if ( z != 0 && z > lastZ) {
            // Consider translation by z only if it has not already been considered
            lastZ = z;
            vector<ActiveTrits> translatedList;
            for (unsigned j = i; j < unitList.size(); j++) {
                ActiveTrits trits = unitList[j]; 
                trits.z = trits.z - z;
                translatedList.push_back(trits); 
            }
            for (unsigned int j = 0; j < i; j++) {
                ActiveTrits trits = unitList[j]; 
                trits.z = trits.z - z + SLICES; 
                translatedList.push_back(trits); 
            }
            // Compare the lists  
            unsigned int idx_cmp; 
            for (idx_cmp = 0; idx_cmp < unitList.size(); idx_cmp++) {
                unsigned int cmp = compare(translatedList[idx_cmp], unitList[idx_cmp]);
                if (cmp == 1)
                    // there is a translated variant smaller than the
                    // original
					return false;
                if (cmp == 2)
                    // original list is smaller than this translated
                    // variant, but still need to check other variants
					break; 
            }
            if (idx_cmp == unitList.size()) {
                // if the two list are identical, then the list is
                // z-periodic, No need to check the other translated 
                // variants.
				break;
            }
        }
    }
    return true; 
}

unsigned int ActiveTritsSet::compare(const ActiveTrits& first, 
                                     const ActiveTrits& second) const
{
    if (first.z < second.z) {
        return 1; 
    } else if (first.z == second.z) {
        if (first.x < second.x) {
            return 1;  
        } else if (first.x == second.x) {
            if (first.isYiActive < second.isYiActive) 
                return 1; 
            else if (first.isYiActive == second.isYiActive) 
                return 0; 
        } 
    }
    return 2; 
}

ActiveTrailCoreCache::ActiveTrailCoreCache()
{
    kernel = true; 
    wA.push(0); 
    wB.push(0);
}
    
MixedTrailCoreCache::MixedTrailCoreCache(const BareState& bareState)
{
    kernel = false;
    activeA.set(bareState.stateA);
    activeB.set(bareState.stateB);
    wA.push(bareState.wA); 
    wB.push(bareState.wB);
    nrComponents = bareState.nrComponents;
    outKernelComponentColumns = bareState.outKernelComponentsColumns;
}

void ActiveTrailCoreCache::pushDummy()
{
    wA.push(wA.top());
	wB.push(wB.top());
	dummy = true;
}

void ActiveTrailCoreCache::push(const ActiveTrits& activeTrits)
{
    dummy = false; 

    unsigned int new_wA = wA.top(); 
    unsigned int new_wB = wB.top(); 

    for (unsigned int y = 0; y < ROWS; y++) {
        if ((activeTrits.isYiActive >> y) & 0x1) {
            TritPosition t(activeTrits.x, y, activeTrits.z);
            if (!activeB.isTheTritInAnActiveTryte(t))
                new_wB += 2; 
            activeB.activateTrit(t); 
            t.invSRSL(); 
            if (!activeA.isTheTritInAnActiveTryte(t))
                new_wA += 2;
            activeA.activateTrit(t); 
        }
    }
    wA.push(new_wA); 
    wB.push(new_wB); 

    if (activeTrits.isYiActive == 0x3 || activeTrits.isYiActive == 0x5 || activeTrits.isYiActive == 0x6) {
        inKernelColumns.push_back(TroikaColumns(activeTrits.x, activeTrits.z)); 
    }

}

void ActiveTrailCoreCache::pop(const ActiveTrits& activeTrits)
{
    wA.pop(); 
    wB.pop(); 
    if (dummy) {
        dummy = false;
        return; 
    }

    for (unsigned int y = 0; y < ROWS; y++) {
        if ((activeTrits.isYiActive >> y) & 0x1) {
            TritPosition t(activeTrits.x, y, activeTrits.z); 
            activeB.deactivateTrit(t); 
            t.invSRSL(); 
            activeA.deactivateTrit(t); 
        }
    }
    if (activeTrits.isYiActive == 0x3 || activeTrits.isYiActive == 0x5 || activeTrits.isYiActive == 0x6) {
        inKernelColumns.pop_back(); 
    }
}
    
void TwoRoundTrailCore::save(ostream& fout)
{
    TroikaState stateB; 
    TroikaState stateA; 
    for (statesB.first(); !statesB.isEnd(); ++statesB) {
        stateB = *statesB; 
        stateA.setInvL(stateB); 
        TrailCore trail(stateA, stateB, wA, wB); 
        trail.save(fout);  
    }
}

void TwoRoundTrailCore::set(const vector<ActiveTrits>& unitList,
                            const ActiveTrailCoreCache& cache, 
                            const TwoRoundTrailCoreCostFunction& costF, 
                            unsigned int aMaxCost)
{
    (void) unitList; 
    (void) costF; 
    wA = cache.wA.top(); 
    wB = cache.wB.top();
    statesB = TroikaStateIterator(cache);
    activeA = cache.activeA;
    activeB = cache.activeB; 

}

void TwoRoundTrailCore::set(const vector<ActiveTrits>& unitList,
                            const MixedTrailCoreCache& cache, 
                            const TwoRoundTrailCoreCostFunction& costF, 
                            unsigned int aMaxCost)
{
    (void) unitList; 
    (void) costF; 
    wA = cache.wA.top(); 
    wB = cache.wB.top();
    statesB = TroikaStateIterator(cache);
    activeA = cache.activeA;
    activeB = cache.activeB; 
}

unsigned int TwoRoundTrailCoreCostFunction::getCost(
                               const vector<ActiveTrits>& unitList,
                               const ActiveTrailCoreCache& cache) const
{
    (void) unitList; 
    return alpha*cache.wA.top() + beta*cache.wB.top();
}

void traverse_K_TrailCoresTree(unsigned int aMaxCost,
                               unsigned int alpha,
                               unsigned int beta, 
                               ostream &fout)
{
    TwoRoundTrailCoreCostFunction costFRun(alpha, beta);  
    ActiveTrailCoreCache activeTritsCache;  
    ActiveTritsSet activeTritsSet;  

    K_TrailCore_Iterator iteratorTrits(activeTritsSet, activeTritsCache, costFRun, aMaxCost, true); 
    for (; !iteratorTrits.isEnd(); ++iteratorTrits) {
        TwoRoundTrailCore trails = *iteratorTrits; 
        trails.save(fout);
    }
}
