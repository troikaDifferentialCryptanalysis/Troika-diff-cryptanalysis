/* backwardInKernelExtension.cpp */
#include "backwardInKernelExtension.h"     
#include "sbox.h"
#include "state.h"

bool InKernelTryteColumns::isEmpty()
{
    return empty; 
}

void InKernelTryteColumns::setFirstValue(const BackwardInKernelExtensionCache& cache)
{
    (void) cache;
    if (empty)
        return; 
    index = 0; 
}

bool InKernelTryteColumns::setNextValue(const BackwardInKernelExtensionCache& cache)
{
    (void) cache;
    if (empty || index >= possibleValues.size() - 1) 
        return false; 
    
    index ++; 
    return true;
}

Weight InKernelTryteColumns::getWeight() const
{
    assert(empty == false); 
    return possibleValues[index].weight; 
}

int InKernelTryteColumns::getHammingWeight() const
{
    assert(empty == false);
    return possibleValues[index].hammingWeight; 
}

bool sortbysecdesc(const pair<int, int> &a,
                   const pair<int, int> &b)
{
       return a.second > b.second;
}

void InKernelTryteColumns::setValues(const TroikaState& state)
{
    tryteAfterST.clear();
    index = 0;

    for (int y = 0; y < ROWS; y++)
        tryteAfterST.push_back(pair<int, int>(y, state.getTryte(xTryte, y, z).value)); 
    sort(tryteAfterST.begin(), tryteAfterST.end(), sortbysecdesc); 
    possibleValues = sbox.inKernelTryteColumnBeforeST[tryteAfterST[0].second][tryteAfterST[1].second][tryteAfterST[2].second];

    if (possibleValues.empty())
        empty = true;
    else 
        empty = false;
}

void InKernelTryteColumns::setTryteColumn(TroikaState &state) const
{
    for (int i = 0; i < 3; i++)
        state.setTryte(xTryte, tryteAfterST[i].first,  z, possibleValues[index].trytes[i]); 
}

ostream & operator << (ostream& fout, const InKernelTryteColumns& col)
{
    fout << "xTryte: " << col.xTryte << " zTryte: " << col.z;
    if (col.empty == false) {
        TroikaState state; 
        col.setTryteColumn(state);
        fout << endl << "TryteColumn " << col.index << " sur "<< col.possibleValues.size() - 1 << " de cout " << 2 *col.getHammingWeight() + (long double)col.getWeight() << endl;
        fout << state; 
    }
    return fout; 
}

BackwardInKernelExtensionPreparation::BackwardInKernelExtensionPreparation
           (long double aMaxWeightExtension, const ActiveState& activeC): 
            maxWeightExtension(aMaxWeightExtension)
{
    int nrActiveTrytesC = 0; 
    map <TrytePosition, int> possibleTrytesAtA;
    TritPosition tA; 

    for (int z = 0; z < SLICES; z++) {
        for (int xTryte = 0; xTryte < 3; xTryte++) {
            int nrActiveTrytesInTheColumn = 0; 
            for (int y = 0; y < ROWS; y++) {
                if (activeC.isTryteActive(xTryte, y, z)) {
                    nrActiveTrytesInTheColumn += 1;
                    for (int tritIndex = 0; tritIndex < 3; ++tritIndex) {
                        tA.set(3 * xTryte + tritIndex, y, z ); 
                        tA.invSRSL(); 
                        possibleTrytesAtA[TrytePosition(tA)] += 1; 
                    }
                }    
            }
            if (nrActiveTrytesInTheColumn == 1) {
                possible = false;
                possibleTrytesAtA.clear();
                return; 
            } else if (nrActiveTrytesInTheColumn > 1)
                tryteColumnsAtB.push_back(InKernelTryteColumns(xTryte, z));
            nrActiveTrytesC += nrActiveTrytesInTheColumn; 
        }
    }
            
    toSubtractFromTheCost = 2 *( 3 * nrActiveTrytesC - possibleTrytesAtA.size());              
    if ( ( 4 * nrActiveTrytesC - toSubtractFromTheCost ) > maxWeightExtension) 
        possible = false; 
     
    else 
        possible = true;
}

void BackwardInKernelExtensionPreparation::prepareExtension(const TroikaState &stateC)
{ 
    vector<InKernelTryteColumns>::iterator it; 
    for (it = tryteColumnsAtB.begin(); it != tryteColumnsAtB.end(); it++) {
        (*it).setValues(stateC);
        if ((*it).isEmpty()) {
            possible = false; 
            return; 
        }
    }
   
    // initialize to add to the cost;
    toAddToTheCost.clear(); 
    toAddToTheCost.push_back(0); 
    Weight toAdd = 0; 
    for (int i = tryteColumnsAtB.size() - 1; i >= 0; i--) {
        toAdd += tryteColumnsAtB[i].getWeight()
                 + 2 * tryteColumnsAtB[i].getHammingWeight();
        toAddToTheCost.insert(toAddToTheCost.begin(), toAdd);
    }
    possible = true;
}

bool BackwardInKernelExtensionPreparation::couldBeExtended() const
{
    return possible; 
}

vector<InKernelTryteColumns>& BackwardInKernelExtensionPreparation::getPartsList() 
{
    return tryteColumnsAtB;
}

ostream & operator << (ostream& fout, const BackwardInKernelExtensionPreparation& prep)
{
    fout << "possible: " << prep.possible << endl; 
    if (prep.possible) {
        fout << "toSubtractFromTheCost: " << prep.toSubtractFromTheCost << endl; 
        fout << "nr columns " << prep.tryteColumnsAtB.size() << endl; 
        for (unsigned int i = 0; i < prep.tryteColumnsAtB.size(); i++) {
            fout << prep.tryteColumnsAtB[i] << endl; 
        }
        fout << "to add to the cost" << endl; 
        for (unsigned int i = 0; i < prep.toAddToTheCost.size(); i++) 
        {
            fout << prep.toAddToTheCost[i] << " ";
        }
        fout << endl; 
    }
    return fout; 
}
 

BackwardInKernelExtensionCache::BackwardInKernelExtensionCache
 (BackwardInKernelExtensionPreparation &prep, const TroikaState& stateC)
:wBC(0), hammingWeightA(0)
{
    prep.prepareExtension(stateC); 
}
 
void BackwardInKernelExtensionCache::push(const InKernelTryteColumns& tryteColumn)
{
    wBC            += tryteColumn.getWeight();
    hammingWeightA += tryteColumn.getHammingWeight();
}

void BackwardInKernelExtensionCache::pop(const InKernelTryteColumns& tryteColumn)
{
    wBC            -= tryteColumn.getWeight();
    hammingWeightA -= tryteColumn.getHammingWeight();
}

CostFunctionBackwardInKernelExtension::CostFunctionBackwardInKernelExtension
(BackwardInKernelExtensionPreparation &prep, const TroikaState& stateC) 
: toAddToTheCost(prep.toAddToTheCost), toSubtractFromTheCost(prep.toSubtractFromTheCost), 
    maxWeightExtension(prep.maxWeightExtension)
{ }

bool CostFunctionBackwardInKernelExtension::tooHighCost(
        const BackwardInKernelExtensionCache& cache, int indCurPart) const
{
    long double cost = 2 * cache.hammingWeightA
                      + (long double) cache.wBC 
                      + (long double)toAddToTheCost[indCurPart + 1] 
                       - toSubtractFromTheCost; 
    if (cost > maxWeightExtension)
        return true; 
    else 
        return false;
}

void BackwardInKernelExtension::set(
                 const vector<InKernelTryteColumns> tryteColumnsAtB,
                 const BackwardInKernelExtensionCache &cache, 
                 const CostFunctionBackwardInKernelExtension &costF)
{
    stateB.setEmptyState(); 
    for (unsigned int i = 0; i < tryteColumnsAtB.size(); i++) 
        tryteColumnsAtB[i].setTryteColumn(stateB); 
    stateA.setInvSRSL(stateB);
    if (costF.toSubtractFromTheCost == 0)
        wMinRevA = 2 * cache.hammingWeightA; 
    else  
        wMinRevA = 2 * stateA.getNrActiveTrytes();
    wBC = cache.wBC; 
}

template<>                     
bool BackwardInKernelExtensionIterator::next()
{
    if (toChild()) {
        if (costF.tooHighCost(cache, indCurPart))
            toParent(); 
        else
            return true;
    }
    do {
        if (toSibling())
            if (!costF.tooHighCost(cache, indCurPart))
                return true; 
        if (!toParent())  
            return false; 
    } while (true); 
}
