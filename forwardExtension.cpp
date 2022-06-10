/* forwardExtension.cpp */ 
#include "forwardExtension.h"
#include "state.h"

TritInfo::TritInfo(unsigned int x, unsigned int y, unsigned int z):
    mustCalculateTheCostOfTheSlice(false),
    mustCalculateTheCostOfTheNextSlice(false), 
    value(0)
{
    posAtSRSLC.set(x, y, z);
    posAtC = posAtSRSLC.getInvSRSL(); 

}

void TritInfo::setFirstValue(const ForwardExtensionCache& cache)
{
    value = 0; 
    // There is always an active value that is valid
    while (!isAValidTritValue(cache.constraintAtC, cache.nrNonActiveTritsAtC))
        value++;
}

bool TritInfo::setNextValue(const ForwardExtensionCache& cache)
{
    unsigned int lastValue = value; 
    do {
        value++; 
    } while ( (value < 3) && ! isAValidTritValue(cache.constraintAtC, cache.nrNonActiveTritsAtC) );
    
    if (value == 3) {
        value = lastValue; 
        return false; 
    }
    return true; 
}

bool TritInfo::isAValidTritValue(
                const map<TritPosition, ConstraintForTritValue>& constraintAtC,
                const map<TrytePosition, unsigned int> &nrNonActiveTritsAtC) const
{

    map<TritPosition, ConstraintForTritValue>::const_iterator constraint;
    map<TrytePosition, unsigned int>::const_iterator nrNonActiveTrits;
    constraint = constraintAtC.find(TritPosition(posAtC));
    nrNonActiveTrits = nrNonActiveTritsAtC.find(TrytePosition(posAtC));
    
    if (value == 0) {
        if (nrNonActiveTrits->second == 2)
            return false;
        if (constraint->second != noConstraint)
            return false; 
    }
    if (value == 1 && constraint->second == mustBe2)
        return false; 
    if (value == 2 && constraint->second == mustBe1)
        return false; 
    return true; 
}

ostream & operator << (ostream &fout, const TritInfo &aTritInfo)
{
    fout << "position at C - SRSL(C) : " << aTritInfo.posAtC << "-" << aTritInfo.posAtSRSLC;
    fout << "value: " << aTritInfo.value; 
    if (aTritInfo.mustCalculateTheCostOfTheSlice)
        fout << " last trit of the silce "; 
    if (aTritInfo.mustCalculateTheCostOfTheNextSlice)
        fout << "next slice passive : " ;
    
    return fout;
}

ForwardExtensionPreparation::ForwardExtensionPreparation(const TroikaState& stateB, 
                                                         long double maxWeightExtension)
:maxWeightExtension(maxWeightExtension)
{
    initPosForSTCompatibility(stateB);
    initPossibleActiveTritsAtSRSLCAndConstraintForActiveTritsAtC(stateB);
    int z; 
    unsigned int zStart; 

    // Choose the slice of the first trit of the vector tritInfoAtSRSLC. 
    for (zStart = 0; zStart < SLICES; zStart++) {
        if (possibleActiveTritsAtSRSLC.isSliceActive(zStart)
                && (!possibleActiveTritsAtSRSLC.isSliceActive((zStart + 1) % SLICES)))
            break; 
    }
    // Add slice by slice the information of the trits of SRSLC
    if (zStart != SLICES) {
        for (unsigned int i = 0; i < SLICES; i++) {
            z = (zStart - i + SLICES) % SLICES; 
            if (possibleActiveTritsAtSRSLC.isSliceActive(z)) {
                addInfoOfTheTritsOfTheSlice(z);  
            }
        }
    }
    else {
        // all the slices are active
        for (z = SLICES - 1; z >= 0; z--) {
            addInfoOfTheTritsOfTheSlice(z);  
            if (z == SLICES - 1)
                tritsInfoAtSRSLC.back().mustCalculateTheCostOfTheSlice = false; 
            if (z == 0)
                tritsInfoAtSRSLC.back().mustCalculateTheCostOfTheNextSlice = true; 
        }
    }
}

void ForwardExtensionPreparation::initPosForSTCompatibility(const TroikaState& stateB)
{
    for (unsigned int z = 0; z < SLICES; z++) {
        for (unsigned xTryte = 0; xTryte < 3; xTryte++) {
            for (unsigned y = 0; y < ROWS; y++) {
                if (stateB.isTryteActive(xTryte, y, z))
                    posForSTCompatibility.push_back(TrytePosition(xTryte, y, z)); 
            }
        }
    }

}

void ForwardExtensionPreparation::
initPossibleActiveTritsAtSRSLCAndConstraintForActiveTritsAtC(const TroikaState &stateB)
{
    TritPosition t; 
    unsigned int tryteValueAtB; 
    ConstraintForTritValue constraint;
    vector<TrytePosition>::const_iterator it; 

    // Use Property 2 of the sbox of Section 3
    for (it = posForSTCompatibility.begin(); 
         it != posForSTCompatibility.end(); it++) {

        tryteValueAtB = (unsigned int)stateB.getTryte(it->x, it->y, it->z);

        // first trit
        t.set(3 * it->x, it->y, it->z);
        constraintAtC[t] = noConstraint;
        t.SRSL(); 
        possibleActiveTritsAtSRSLC.activateTrit(t);

        // second trit 
        if (tryteValueAtB == 9 || tryteValueAtB == 18)
            constraint = cannotBe0; 
        else 
            constraint = noConstraint; 
        t.set(3 * it->x + 1, it->y, it->z);
        constraintAtC[t] = constraint;
        t.SRSL(); 
        possibleActiveTritsAtSRSLC.activateTrit(t);

        // third trit
        if (tryteValueAtB == 1)
            constraint = mustBe1; 
        else if (tryteValueAtB == 2)
            constraint = mustBe2; 
        else 
            constraint = noConstraint; 
        t.set(3 * it->x + 2, it->y, it->z); 
        constraintAtC[t] = constraint;
        t.SRSL(); 
        possibleActiveTritsAtSRSLC.activateTrit(t);

    }
}

void ForwardExtensionPreparation::addInfoOfTheTritsOfTheSlice(unsigned int z)
{ 
    for (unsigned int x = 0; x < COLUMNS; x++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            if (possibleActiveTritsAtSRSLC.isTritActive(x, y, z)) {
                tritsInfoAtSRSLC.push_back(TritInfo(x, y, z));
            }
        }
    }
    tritsInfoAtSRSLC.back().mustCalculateTheCostOfTheSlice = true;
    if (!possibleActiveTritsAtSRSLC.isSliceActive((z - 1 + SLICES) % SLICES))
        tritsInfoAtSRSLC.back().mustCalculateTheCostOfTheNextSlice = true; 
}

ostream & operator << (ostream& fout, const ForwardExtensionPreparation& prep)
{

    fout << "posForSTCompatibility: " << endl; 
    for (unsigned int i = 0; i < prep.posForSTCompatibility.size(); i++) {
        cout << prep.posForSTCompatibility[i] << endl; 
    }
    fout << "constraints at C: " << endl;
    fout << "possible active trits at SRSLC: " << endl; 
    fout << prep.possibleActiveTritsAtSRSLC << endl; 
    map<TritPosition, ConstraintForTritValue>::const_iterator it; 
    for (it = prep.constraintAtC.begin(); it != prep.constraintAtC.end(); ++it) {
        fout << it->first << " " << it->second << endl; 
    }
    fout << "tritsInfoAtSRSLC: " << endl; 
    for (unsigned int i = 0; i < prep.tritsInfoAtSRSLC.size(); i++) {
        cout << prep.tritsInfoAtSRSLC[i] << endl; 
    }
    return fout; 
}

ForwardExtensionCache::ForwardExtensionCache(ForwardExtensionPreparation& prep, 
    const TroikaState& stateB):
constraintAtC(prep.constraintAtC)
{
    (void) stateB;
    stackNrActiveTrytesD.push(0); 
    for (unsigned int x = 0; x < COLUMNS; x++) {
        for (unsigned int z = 0; z < SLICES; z++) {
            columnParityOfSRSLC[x][z] = 0; 
            for (unsigned int y = 0; y < ROWS; y++) {
                tritsAtSRSLC[x][y][z] = 0; 
                tritsAtD[x][y][z] = 0; 
            }
        }
    } 

}

void ForwardExtensionCache::push(const TritInfo &trit)
{
    TritPosition t = trit.posAtSRSLC;
    tritsAtSRSLC[t.x][t.y][t.z] = trit.value;
    
    if (trit.value == 0)
        nrNonActiveTritsAtC[TrytePosition(trit.posAtC)] += 1; 

    if (trit.mustCalculateTheCostOfTheSlice) {
        setTheColumnParityOfTheSliceAtSRSLC(t.z);
        applyToTheSliceAddColumnParity(t.z);
        stackNrActiveTrytesD.push(stackNrActiveTrytesD.top()
                                     + getNrActiveTrytesOfTheSliceAtD(t.z));
    }
    if (trit.mustCalculateTheCostOfTheNextSlice) {
        unsigned int z = (t.z - 1 + SLICES) % SLICES ;
        applyToTheSliceAddColumnParity(z); 
        stackNrActiveTrytesD.push(stackNrActiveTrytesD.top()
                                     + getNrActiveTrytesOfTheSliceAtD(z));
    } 

}

void ForwardExtensionCache::pop(const TritInfo &trit)
{
    TritPosition t = trit.posAtSRSLC;
    tritsAtSRSLC[t.x][t.y][t.z] = 0;

    if (trit.value == 0)
        nrNonActiveTritsAtC[TrytePosition(trit.posAtC)] -= 1; 

    if (trit.mustCalculateTheCostOfTheNextSlice) {
        stackNrActiveTrytesD.pop();
    }
    if (trit.mustCalculateTheCostOfTheSlice) {
        stackNrActiveTrytesD.pop();
    } 
}
     
void ForwardExtensionCache::setTheColumnParityOfTheSliceAtSRSLC(unsigned int z)
{
    unsigned int columnParity;
    for (unsigned int x = 0; x < COLUMNS; x++) {
        columnParity = 0; 
        for (unsigned int y = 0; y < ROWS; y++) {
            columnParity += tritsAtSRSLC[x][y][z]; 
        }
        columnParityOfSRSLC[x][z] = columnParity % 3; 
    }
}

void ForwardExtensionCache::applyToTheSliceAddColumnParity(unsigned int z)
{
    for (unsigned int x = 0; x < COLUMNS; x++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            
            tritsAtD[x][y][z] 
                 = tritsAtSRSLC[x][y][z] 
                 + (columnParityOfSRSLC[(x - 1 + COLUMNS) % COLUMNS][z]
                 + columnParityOfSRSLC[(x + 1) % COLUMNS][( z + 1) % SLICES] ); 
            tritsAtD[x][y][z] %= 3; 
        }
    }

}

unsigned int ForwardExtensionCache::getNrActiveTrytesOfTheSliceAtD(
                                                        unsigned int z) const
{
    unsigned int nrActiveTrytesAtD = 0; 
    for (unsigned int x_tryte = 0; x_tryte < 3; x_tryte++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            if (   tritsAtD[3 * x_tryte][y][z]
                || tritsAtD[3 * x_tryte + 1][y][z]
                || tritsAtD[3 * x_tryte + 2][y][z] )
                nrActiveTrytesAtD += 1; 
        }
    }
    return nrActiveTrytesAtD; 
}

ostream & operator << (ostream& fout, const ForwardExtensionCache& cache)
{
    for (map<TritPosition, ConstraintForTritValue>::const_iterator it = cache.constraintAtC.begin(); 
            it != cache.constraintAtC.end(); it++) {
        cout << it->first << " "; 
        if (it->second == noConstraint)
            cout << "noConstraint" << endl; 
        if (it->second == mustBe1)
            cout << "mustBe1" << endl; 
        if (it->second == mustBe2)
            cout << "mustBe2" << endl; 
        if (it->second == cannotBe0)
            cout << "cannotBe0" << endl; 
    }
    for (map<TrytePosition, unsigned int>::const_iterator it = cache.nrNonActiveTritsAtC.begin(); 
            it != cache.nrNonActiveTritsAtC.end(); it++) {
        cout << it->first << " " << it->second << endl; 
    }
    TroikaState stateSRSLC; 
    TroikaState stateD; 
    for (unsigned int x = 0; x < COLUMNS; x++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            for (unsigned int z = 0; z < SLICES; z++) {
                stateD.setTritValue(cache.tritsAtD[x][y][z], x, y, z);
                stateSRSLC.setTritValue(cache.tritsAtSRSLC[x][y][z], x, y, z);
            }
        }
    }
    fout << "nrActiveTrytesAtD: " << cache.stackNrActiveTrytesD.top() << endl;

    for (unsigned int z = 0; z < SLICES; z++) {
        bool nonzeroPartity = false; 
        for (unsigned int x = 0; x < COLUMNS; x++) {
            if (cache.columnParityOfSRSLC[x][z]) 
                nonzeroPartity = true; 
        }
        if (nonzeroPartity) {
            cout <<  setw(4) << "z: "  << z << " "; 
            for (unsigned int x = 0; x < COLUMNS; x++) 
                fout << cache.columnParityOfSRSLC[x][z] << " ";
            fout << endl; 
        }
    }
    fout << "state SRSLC : " << endl << stateSRSLC << endl; 
    fout << "active state D: " << endl << stateD << endl; 
    fout << cache.activeD << endl;
    return fout; 
}

CostFunctionForwardExtension::CostFunctionForwardExtension(ForwardExtensionPreparation &prep, 
                                                           const TroikaState& stateC)
{
    (void) stateC; 
    maxWeightExtension = prep.maxWeightExtension; 
    nrActiveTrytesC = prep.posForSTCompatibility.size();
}

bool CostFunctionForwardExtension::tooHighCost(const ForwardExtensionCache& cache,
                                               int indCurPart) const
{
    if ( 2 * (cache.stackNrActiveTrytesD.top() + nrActiveTrytesC) > maxWeightExtension ) 
        return true; 
    return false; 
}

ForwardExtension::ForwardExtension(const ForwardExtensionPreparation& prep, const TroikaState &stateB)
:stateB(stateB), posForSTCompatibility(prep.posForSTCompatibility){}


ForwardExtension::ForwardExtension(const TroikaState& aStateB, const vector<TrytePosition> posForSTCompatibility)
: stateB(aStateB), posForSTCompatibility(posForSTCompatibility){}

bool ForwardExtension::isValidAndBelowWeight(long double maxWeightExtension) const
{
    if ( valid && ((long double)getWeight() <= maxWeightExtension) ) 
        return true; 
    return false; 
}

void ForwardExtension::setStateCAndDFromStateD(const TroikaState& aStateD)
{
    stateD = aStateD; 
    stateC.setInvL(aStateD);
    valid = sbox.areSTCompatible(stateB, stateC, posForSTCompatibility, wBC);
    wMinDirD =  2 * stateD.getNrActiveTrytes();
}  


 
void ForwardExtension::set(const vector<TritInfo>& tritsInfoAtSRSLC, 
                           const ForwardExtensionCache& cache, 
                           const CostFunctionForwardExtension& costF)
{
    // initialize stateC and stateD 
    stateD.setEmptyState();
    vector<TritInfo>::const_iterator it; 
    for (it = tritsInfoAtSRSLC.begin(); it != tritsInfoAtSRSLC.end(); it++) {
        stateD.setTritValue(it->value, it->posAtSRSLC.x, it->posAtSRSLC.y, it->posAtSRSLC.z); 
    }
    stateD.addColumnParity(); 
    stateC.setInvL(stateD); 
    wMinDirD = 2 * cache.stackNrActiveTrytesD.top();
    valid = sbox.areSTCompatible(stateB, stateC, posForSTCompatibility, wBC);
}

Weight ForwardExtension::getNrActiveTrytesC() const
{
    return posForSTCompatibility.size();
}
 
template<>
bool ForwardExtensionIterator::next()
{
    if (toChild()) {
        if (false == costF.tooHighCost(cache, indCurPart))
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
