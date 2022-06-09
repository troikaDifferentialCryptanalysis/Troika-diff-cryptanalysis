/* backwardExtension.cpp */
#include "backwardExtension.h"
#include "state.h"

TryteInfo::TryteInfo(const TrytePosition &aPosition, const Tryte & tryteAtC)
{
    Sbox sbox; 
    position = aPosition; 
    possibleValues = sbox.inputDiff[(unsigned int) tryteAtC];
    index = 0; 
    mustCalculateTheCostOfTheSlice = false; 
    mustCalculateTheCostOfTheNextSlice = false; 
}

void TryteInfo::setFirstValue(const BackwardExtensionCache& cache)
{ 
    (void) cache; 
    index = 0; 
}

bool TryteInfo::setNextValue(const BackwardExtensionCache& cache)
{
    if (index >= possibleValues.size() - 1) 
        return false; 
    index ++; 
    return true;
}

Tryte TryteInfo::getValue() const
{
    assert(index < possibleValues.size());
    return possibleValues[index];
}

Weight TryteInfo::getWeight() const
{
    assert(index < possibleValues.size());
    return possibleValues[index].weight;
}

ostream & operator << (ostream &fout, const TryteInfo &aTryteInfo)
{
    fout << aTryteInfo.position; 
    fout << "index  " << aTryteInfo.index << " sur " << aTryteInfo.possibleValues.size() - 1; 
    if (aTryteInfo.mustCalculateTheCostOfTheSlice)
        fout << " endSlice  "; 
    if (aTryteInfo.mustCalculateTheCostOfTheNextSlice)
        fout << " next slice passive "; 
    return fout;
}

BackwardExtensionPreparation::BackwardExtensionPreparation(
                                            const TroikaState &stateC, 
                                            long double maxWeightExtension)
:maxWeightExtension(maxWeightExtension)
{
    // Initialize the vector trytesInfoAtB
    bool sliceActive[SLICES];
    bool allSlicesActive = true; 
    for (unsigned int z = 0; z < SLICES; z++) {
        if (stateC.isSliceActive(z))
            sliceActive[z] = true; 
        else {
            sliceActive[z] = false; 
            allSlicesActive = false;  
        }
    }
    if (allSlicesActive == true) {
        for (int z = SLICES - 1; z >= 0; z--) {
            addInfoOfTheTrytesOfTheSliceAtB(z, stateC); 
            if (z == SLICES - 1)
                trytesInfoAtB.back().mustCalculateTheCostOfTheSlice = false; 
            if (z == 0)
                trytesInfoAtB.back().mustCalculateTheCostOfTheNextSlice = true; 
        }
        return;
    } 

    vector<pair<unsigned int, unsigned int> > firstZConsecutiveActiveSlices;
    for (unsigned int z = 0; z < SLICES; z++) {
        if (sliceActive[z] && ! sliceActive[(z + 1) % SLICES]) {
            unsigned int nrActiveTrytes = stateC.getNrActiveTrytesOfSlice(z);
            pair<unsigned int, unsigned int> nrActiveTrytesOfSlices(nrActiveTrytes, z);
            firstZConsecutiveActiveSlices.push_back(nrActiveTrytesOfSlices);
        }
    }
    sort(firstZConsecutiveActiveSlices.begin(), firstZConsecutiveActiveSlices.end());

    vector<pair<unsigned int, unsigned int> >::const_iterator it; 
    for (it = firstZConsecutiveActiveSlices.begin(); it != firstZConsecutiveActiveSlices.end(); ++it) {
        int z = it->second;
        do {
            addInfoOfTheTrytesOfTheSliceAtB(z, stateC);  
            z = (z - 1 + SLICES) % SLICES; 
        } while (sliceActive[z] == true);
    }
}

vector<TryteInfo>& BackwardExtensionPreparation::getPartsList()
{
    return trytesInfoAtB;
}

void BackwardExtensionPreparation::addInfoOfTheTrytesOfTheSliceAtB(unsigned int z,
                                                                   const TroikaState &stateC)
{
    
    for (unsigned int xTryte = 0; xTryte < 3; xTryte++) {
        for (unsigned y = 0; y < ROWS; y++) {
            if (stateC.isTryteActive(xTryte, y, z)) {
                Tryte tryteValueAtC = stateC.getTryte(xTryte, y, z); 
                TrytePosition pos(xTryte, y, z); 
                trytesInfoAtB.push_back(TryteInfo(pos, tryteValueAtC));
            } 
        }
    }
    trytesInfoAtB.back().mustCalculateTheCostOfTheSlice = true;
    if (!stateC.isSliceActive((z - 1 + SLICES) % SLICES))
        trytesInfoAtB.back().mustCalculateTheCostOfTheNextSlice = true; 
    
}

ostream & operator << (ostream &fout, const BackwardExtensionPreparation& prep)
{
    fout << prep.trytesInfoAtB.size() << " nr active trytes " << endl; 
    for (unsigned int i = 0; i < prep.trytesInfoAtB.size(); i++)
        fout << prep.trytesInfoAtB[i] << endl;
    return fout; 
}

BackwardExtensionCache::BackwardExtensionCache(const BackwardExtensionPreparation& prep, 
                                               const TroikaState& stateC)
{
    (void) prep;
    (void) stateC;
    stackWeightBC.push(0); 
    nrActiveTrytesA = 0; 
    for (unsigned int x = 0; x < COLUMNS; x++) {
        for (unsigned int z = 0; z < SLICES; z++) {
            columnParityOfB[x][z] = 0; 
            for (unsigned int y = 0; y < ROWS; y++) {
                B[x][y][z] = 0; 
                invAddColumnParityOfB[x][y][z] = 0; 
            }
        }
    }    
}

void BackwardExtensionCache::push(const TryteInfo& tryteAtB)
{ 
    stackWeightBC.push(stackWeightBC.top() + tryteAtB.getWeight());
    initTryteOfStateB(tryteAtB);
    if (tryteAtB.mustCalculateTheCostOfTheSlice) {
        unsigned int z = tryteAtB.position.z; 
        setTheColumnParityOfTheSlice(z);
        applyToTheSliceInvAddColumnParity(z);
        addTrytesAtAFromASliceOfInvAddColumnParityOfB(z); 
    }
    if (tryteAtB.mustCalculateTheCostOfTheNextSlice) {
        unsigned int z = (tryteAtB.position.z - 1 + SLICES) % SLICES ;
        applyToTheSliceInvAddColumnParity(z); 
        addTrytesAtAFromASliceOfInvAddColumnParityOfB(z);
    }
}
        
void BackwardExtensionCache::pop(const TryteInfo& tryteAtB)
{
    // It is not necessary to update the tryte value. It is going to be updated 
    // by the next call to push.
    if (tryteAtB.mustCalculateTheCostOfTheNextSlice) {
        unsigned int z = (tryteAtB.position.z - 1 + SLICES) % SLICES ; 
        removeTrytesAtAFromASliceOfInvAddColumnParityOfB(z);
    }
    if (tryteAtB.mustCalculateTheCostOfTheSlice) 
        removeTrytesAtAFromASliceOfInvAddColumnParityOfB(tryteAtB.position.z); 
    stackWeightBC.pop();
}
 
ostream & operator << (ostream &fout, const BackwardExtensionCache& cache)
{
    TroikaState stateB; 
    TroikaState stateInvAddColParityB; 
    for (unsigned int x = 0; x < COLUMNS; x++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            for (unsigned int z = 0; z < SLICES; z++) {
                stateB.setTritValue(cache.B[x][y][z], x, y, z);
                stateInvAddColParityB.setTritValue(cache.invAddColumnParityOfB[x][y][z], x, y, z);
            }
        }
    }
    fout << "wBC: " << cache.stackWeightBC.top() << endl;

    for (unsigned int z = 0; z < SLICES; z++) {
        bool nonzeroPartity = false; 
        for (unsigned int x = 0; x < COLUMNS; x++) {
            if (cache.columnParityOfB[x][z]) 
                nonzeroPartity = true; 
        }
        if (nonzeroPartity) {
            cout << setw(4) << "z: "  << z << " "; 
            for (unsigned int x = 0; x < COLUMNS; x++) 
                fout << cache.columnParityOfB[x][z] << " ";
            fout << endl; 
        }
    }
    fout << "state B : " << endl << stateB << endl; 
    fout << "state invAddColumnParityOfB: " << endl << stateInvAddColParityB << endl; 
    fout << cache.nrActiveTrytesA << " active trytes at A: " << endl; 
    fout << cache.activeA << endl;
    return fout; 
}
 
void BackwardExtensionCache::initTryteOfStateB(const TryteInfo& tryte)
{
    for (int i = 0; i < 3; i++) {
        B[3 * tryte.position.x + i][tryte.position.y][tryte.position.z]
        = tryte.getValue().trit(i);
    }
}
     
void BackwardExtensionCache::applyToTheSliceInvAddColumnParity(unsigned int z)
{
    for (unsigned int x = 0; x < COLUMNS; x++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            int colParity_1 = columnParityOfB[(x - 1 + COLUMNS) % COLUMNS][z];
            int colParity_2 = columnParityOfB[(x + 1) % COLUMNS][ (z + 1) % SLICES ];
            invAddColumnParityOfB[x][y][z] = B[x][y][z] + 2 * colParity_1 + 2 * colParity_2; 
            invAddColumnParityOfB[x][y][z] %= 3; 
        }
    }
}

void BackwardExtensionCache::setTheColumnParityOfTheSlice(unsigned int z)
{
    unsigned int columnParity; 
    for (unsigned int x = 0; x < COLUMNS; x++) {
        columnParity = 0; 
        for (unsigned int y = 0; y < ROWS; y++) {
            columnParity += B[x][y][z];
        }
        columnParityOfB[x][z] = columnParity % 3; 
    }
}

void BackwardExtensionCache::addTrytesAtAFromASliceOfInvAddColumnParityOfB(unsigned int z) 
{
    TritPosition tA; 
    for (unsigned int x = 0; x < COLUMNS; x++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            if (invAddColumnParityOfB[x][y][z] != 0) {
                tA.set(x, y, z); 
                tA.invSRSL();
                if (activeA.isTheTritInAnActiveTryte(tA) == false)
                    nrActiveTrytesA++;
                activeA.activateTrit(tA); 
            }
        }
    }
}

void BackwardExtensionCache::removeTrytesAtAFromASliceOfInvAddColumnParityOfB(unsigned int z)
{
    TritPosition tA; 
    for (unsigned int x = 0; x < COLUMNS; x++) { 
        for (unsigned int y = 0; y < ROWS; y++) { 
            if (invAddColumnParityOfB[x][y][z] != 0) {
                tA.set(x, y, z); 
                tA.invSRSL();
                activeA.deactivateTrit(tA); 
                if (activeA.isTheTritInAnActiveTryte(tA) == false)
                    nrActiveTrytesA --; 
            }
        }
    }
}
 
CostFunctionBackwardExtension::CostFunctionBackwardExtension(
                                      const BackwardExtensionPreparation& prep, 
                                      const TroikaState& stateC):
maxWeightExtension(prep.maxWeightExtension), nrActiveTrytesB(prep.trytesInfoAtB.size()){}

bool CostFunctionBackwardExtension::tooHighCost(
        const BackwardExtensionCache& cache, int indCurPart) const
{
    long double cost = 2 * (nrActiveTrytesB - indCurPart - 1)
                       + (long double) cache.stackWeightBC.top()
                       + 2 *  cache.nrActiveTrytesA;  
    if (cost > maxWeightExtension)
        return true; 
    else 
        return false;
}

void BackwardExtension::set(const vector<TryteInfo> trytesInfoAtB,
         const BackwardExtensionCache& cache, 
         const CostFunctionBackwardExtension& costF)
{
    stateB.setEmptyState();
    for (unsigned int i = 0; i < trytesInfoAtB.size(); i++) {
        stateB.setTryte(trytesInfoAtB[i].position.x,
                          trytesInfoAtB[i].position.y,
                          trytesInfoAtB[i].position.z,
                          trytesInfoAtB[i].getValue()); 
    }    
    stateA.setInvL(stateB);
    wBC = cache.stackWeightBC.top();
    wMinRevA = 2 * cache.nrActiveTrytesA;

}

bool BackwardExtension::isValidAndBelowWeight(long double maxWeightExtension) const
{
    if ( (long double)(wBC + wMinRevA)<= maxWeightExtension ) 
        return true;  
    return false; 
}

    template<>
bool BackwardExtensionIterator::next()
{
    if (toChild()) {
        if (! costF.tooHighCost(cache, indCurPart))
            return true; 
    }

    do {
        while (toSibling()) {
            if (! costF.tooHighCost(cache, indCurPart))
                return true;
        }
        if (!toParent()) 
            return false; 
    } while (true); 
}
