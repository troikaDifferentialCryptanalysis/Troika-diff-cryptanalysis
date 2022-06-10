/** bareStateIterator.cpp */ 
#include "bareStateIterator.h"
#include "state.h"
#include "traversal.h"
#include "troikaStateIterator.h"
#include "types.h"
#include <iterator>
#include <vector>

/** Values of the column of an affected unit before theta. */ 
const array<int, 3> ZERO_PARITY_COLUMN[9] =
{  {0, 0, 0}, {1, 2, 0}, {2, 1, 0}, 
   {1, 0, 2}, {2, 0, 1}, {0, 1, 2}, 
   {0, 2, 1}, {1, 1, 1}, {2, 2, 2}};

Column::Column()
{
    // initialize the column with the smallest column
    dSupraUnit = 0;  
    zSupraUnit = 0; 
    rank       = 0;
    ending     = false;
    indexValue = 0; 
    scalar     = 1; 
    entanglementType = noEntanglement; 
    setCoordinates(); 
}

void Column::incrementRank()
{
    rank++; 
    indexValue = 0; 
    setCoordinates();
}

bool Column::incrementIndexValue()
{
    if (affectedBy() != 0) {
        if (indexValue < 8)
            indexValue ++; 
        else 
            return false; 
    } else { 
        // Only unaffected unit of indexValue == 0 can be concern be an overlap.
        if (entanglementType == noEntanglement && indexValue < 2)
            indexValue ++; 
        else
            return false;
    }
    return true; 
}

bool Column::incrementSupraUnitPosition(const vector<Column>& unitList) 
{

    if (unitList.empty()) {
        // optimization for z-canonicity : the coordinates of the first
        // supra-unit must satisfy dSupraUnit = 0 and 0 ≤ zSupraUnit < 9
        if (zSupraUnit < 8) {
            zSupraUnit++; 
            indexValue = 0; 
        } 
        else 
            return false; 
    } else {
        // optimization for z-canonicity : a supra-unit translated along the 
        // z-axes to have dSupraUnit = 0 cannot come before the first supra-unit 
        // of the list. 
        unsigned int zFirstSupraUnit = unitList[0].zSupraUnit;
        do {
            if (zSupraUnit < SLICES - 1) 
                zSupraUnit++; 
            else if (dSupraUnit < DIAGONAL - 1) {
                dSupraUnit++; 
                zSupraUnit = 0; 
            } else 
                return false; 
        } while (zFirstSupraUnit > zSupraUnitCanonical()); 
        
        if (zFirstSupraUnit == zSupraUnitCanonical())
            indexValue = unitList[0].indexValue; 
        else 
            indexValue = 0;

    }
    setCoordinates(); 
    return true; 
}

void Column::changeScalar()
{
    scalar = (2 * scalar) % 3; 
}


unsigned int Column::affectedBy() const
{
    if (rank % 2 == 1)
        return 0; 
    return scalar; 
}

unsigned int Column::parity() const
{
    if (rank % 2 == 0)
        return 0; 
    return scalar; 
}

unsigned int Column::zSupraUnitCanonical() const
{
    return (5*dSupraUnit + zSupraUnit) % 9; 
} 

void Column::zTranslateSupraUnit(int dz)
{
    unsigned int xSupraUnit = (dSupraUnit + 2 * zSupraUnit) % COLUMNS;
    zSupraUnit = (( (zSupraUnit + dz ) % SLICES + SLICES ) % SLICES);      
    dSupraUnit = (xSupraUnit - 2 * zSupraUnit + 2 * SLICES) % COLUMNS;
}

bool Column::operator < (const Column& other) const
{
    if (dSupraUnit < other.dSupraUnit)
        return true; 
    else if (dSupraUnit == other.dSupraUnit) {
        if (zSupraUnit < other.zSupraUnit) 
            return true; 
        else if (zSupraUnit == other.zSupraUnit) {
            if (rank < other.rank)
                return true; 
            else if (rank == other.rank) {
                if ((ending == false) && (other.ending == true))
                    return true; 
                else if (ending == other.ending) {
                    if (indexValue < other.indexValue)
                        return true; 
                }
            }
        }
    }
    return false; 
}

ostream & operator << (ostream& fout, const Column& col)
{
    fout << "(" << dec << col.x << ", " << col.z << ", ";
    fout << " r = " << col.rank << ", "; 
    if (col.affectedBy())
        fout << "Aff by " << col.affectedBy() << ", "; 
    if (col.parity())
        fout << "Parity " << col.parity() << ", "; 
    if (col.ending)
        fout << "ending, " ; 
    fout << "index " << col.indexValue << ", ";
    fout << "scalar " << col.scalar << ", "; 
    if (col.entanglementType == noEntanglement)
        fout << "noEntanglement" ; 
    if (col.entanglementType == noNeedToIterate)
        fout << "noNeedToIterate" ; 
    if (col.entanglementType == toIterate)
        fout << "toIterate" ; 
    if (col.entanglementType == iterated)
        fout << "iterated" ;
    fout << " )";
    fout << " dSupraUnit: " << col.dSupraUnit; 
    fout << " zSupraUnit: " << col.zSupraUnit; 
    return fout; 
}

void Column::setCoordinates()
{
    unsigned int d;  
    if (rank % 2 == 0) {
        d = (dSupraUnit + 1) % 9; 
        z = (SLICES + zSupraUnit - 1 + rank/2) % SLICES;
    } else {
        d = dSupraUnit; 
        z = (SLICES + zSupraUnit + rank/2 ) % SLICES; 
    }
    x = (d + 2 * z) % COLUMNS; 
}

Column ColumnsSet::getFirstChildUnit(vector<Column>& unitList,
                                     BareStateCache& cache) const
{
 
    if (unitList.empty())
        return Column();
    const Column& parent = unitList.back();
    Column child = parent;

    // ENDING COLUMN --> STARTING COLUMN 
    if (parent.ending) {
        child.rank   = 0;
        child.ending = false;
        child.scalar = 1;

        do {
            if (!child.incrementSupraUnitPosition(unitList))
                throw EndOfSet(); 
        } while (checkColumnOverlapping(unitList, child, cache)); 
         
    }
    // CONTINUE THE SAME SUPRA UNIT
    else {
        child.incrementRank();
        if (parent.rank != 0)
            child.changeScalar();
        if (checkColumnOverlapping(unitList, child, cache)) {
            if (parent.parity() != 0) {
                child.incrementRank(); 
                if (checkColumnOverlapping(unitList, child, cache))
                    throw EndOfSet();
            } else 
                throw EndOfSet();
        }
    }
    return child;
}

void ColumnsSet::iterateUnit(vector<Column>& unitList, Column& current,
                             BareStateCache& cache) const
{
    // STARTING COLUMN 
    if (current.rank == 0) {
        if (!current.incrementIndexValue()) {
            if (current.entanglementType == toIterate) 
                cache.iterateEntanglement(unitList, current);
            else {
                // Move the starting column
                do {
                    if (!current.incrementSupraUnitPosition(unitList))
                        throw EndOfSet();

                } while (checkColumnOverlapping(unitList, current, cache));
            }
        }

    } 
    // UNAFFECTED COLUMN 
    else if (current.affectedBy()== 0) {
        if (!current.incrementIndexValue()) {
            if (current.entanglementType == toIterate)
                cache.iterateEntanglement(unitList, current); 
            else 
                throw EndOfSet();      
        }
    }
    // AFFECTED MIDDLE AND ENDING COLUMN  
    else {
        if (!current.incrementIndexValue()) {
            if (current.ending == false) {
                current.indexValue = 0; 
                current.ending = true; 
                current.changeScalar(); 
            } else {
                if (current.entanglementType == toIterate) {
                    cache.iterateEntanglement(unitList, current); 
                    current.ending = false; 
                    current.changeScalar(); 
                } else {
                    current.ending = false; 
                    current.incrementRank(); 
                    current.changeScalar(); 
                    if (checkColumnOverlapping(unitList, current, cache))
                        throw EndOfSet(); 
                }
            }
        }
    }
}

bool ColumnsSet::isCanonical(const vector<Column>& unitList, 
                             const BareStateCache& cache) const
{
    if (unitList[0].dSupraUnit != 0 || unitList[0].zSupraUnit > 8)
        return false; 
    if (unitList.back().ending == false)
        return true;
    if (cache.startSupraUnit.size() == 0)
        return true; 
    vector<unsigned int>::const_iterator ind = cache.startSupraUnit.begin();
    for (ind++; ind != cache.startSupraUnit.end(); ind++) {

        int dz = unitList[*ind].zSupraUnitCanonical() - unitList[*ind].zSupraUnit;
        // create a translated unit-list and sort its element
        vector<Column> translatedList = unitList; 
        for (vector<Column>::iterator col = translatedList.begin(); col != translatedList.end(); col++)
            (*col).zTranslateSupraUnit(dz); 
        sort(translatedList.begin(), translatedList.end());
        // compare the initial list to the translated list
        for (unsigned int idx_cmp = 0; idx_cmp < unitList.size(); idx_cmp++) {
            if (translatedList[idx_cmp] < unitList[idx_cmp])
                // there is a translated variant smaller than the 
                // original
                return false;
            if (unitList[idx_cmp] < translatedList[idx_cmp])
                // original list is smaller than this translated
                // variant, but still need to check other variant
                break;
        }
    }
    return true; 
}

bool ColumnsSet::checkColumnOverlapping(const vector<Column>& unitList,
                                        Column& current,
                                        const BareStateCache& cache) const
{

    int colOfSameXZ = cache.supraUnitIndexes[current.x][current.z]; // not unsigne int
    if (colOfSameXZ == -1) {
        current.entanglementType = noEntanglement; 
        return false; 
    }

    // There is an overlap. We must verify if it is authorized or not.  
    if (current.affectedBy() != 0) {
        // An affected column cannot overlap another affected column
        if (cache.thetaEffect.isTritActive(current.x, current.z))
            return true;
        // An affected column can only overlap a column with a single trit at y = 0
        if (false == cache.stateB.isTritActive(current.x, 0, current.z))
            return true; 
        else {
            // The overlap is authorized. 
            if (current.rank == 0)
                current.entanglementType = toIterate;
            else if (cache.lastSupraUnitInTheComponentOf(colOfSameXZ)) {
                current.entanglementType = noNeedToIterate; 
            } else 
                current.entanglementType = toIterate;
        }
    }
    else {
        // An unaffected column cannot overlap another unaffected column
        if (cache.parityPlane.isTritActive(current.x, current.z))
            return true;
        // The overlap is authorized.
        if (cache.lastSupraUnitInTheComponentOf(colOfSameXZ)) {
            current.entanglementType = noNeedToIterate; 
        } else 
            current.entanglementType = toIterate;
    }
    return false; 
}

ostream & operator << (ostream& fout, const TritPositionAtAAndB& trit)
{
    fout << trit.pA << " - " << trit.pB; 
    return fout;
}

TritAtAAndB::TritAtAAndB(const Column& col, unsigned int y):
    TritPositionAtAAndB(col, y)
{ 
    // initialize the trit values
    if (col.affectedBy() == 0) {
        vB = (y == col.indexValue ?  1 : 0); 
        vA = vB; 
    } else {
        vA = ZERO_PARITY_COLUMN[col.indexValue][y]; 
        vB = vA + 1; 
    }
    vA  = (vA * col.scalar) % 3; 
    vB  = (vB * col.scalar) % 3;

    dSupraUnit = col.dSupraUnit; 

    // initialize the attribute stable 
    // (see paragraph "importance of the ordering of supra-units", Section 5.4)
    if (y != 0) 
        stable = true; 
    else if (col.entanglementType != noEntanglement)
        stable = true; 
    else if (vA == vB) { // If the unit is unaffected 
        if (dSupraUnit > 0)
            stable = true; 
        else 
            stable = false; 
    } else {
        stable = false; 
    }
}

ostream & operator << (ostream& fout, const TritAtAAndB& trit)
{
    fout << trit.pA << ": " << trit.vA; 
    fout << " - " << trit.pB << ": " << trit.vB;
    if (trit.stable)
        fout << " stable "; 
    else 
        fout << " unstable"; 
    return fout;
}

BareStateCache::BareStateCache():
    nrActiveTrytesA(0), nrActiveTrytesB(0), nrStableTrytesA(0), 
    nrStableTrytesB(0), dummy(false), indexLastSupraUnit(-1), 
    indexLastUnit(-1)
{
    for (unsigned x = 0; x < COLUMNS; x++) {
        for (unsigned int z = 0; z < SLICES; z++) {
            supraUnitIndexes[x][z] = -1; 
        }
    }
    indexLastSupraUnit = -1;
    indexLastUnit = -1; 
}
        
void BareStateCache::pushDummy()
{
    dummy = true; 
}

void BareStateCache::push(const Column& col)
{
    indexLastUnit++; 
    thetaEffect.addTritValue(col.affectedBy(), col.x, col.z); 
    parityPlane.addTritValue(col.parity(), col.x, col.z); 

    if (col.affectedBy() != 0) {
        for (unsigned int y = 0; y < ROWS; y++) 
            pushOrPopTritAtAAndB(true, TritAtAAndB(col, y));
    } else {
        pushOrPopTritAtAAndB(true, TritAtAAndB(col, col.indexValue));
    }

    if (col.rank == 0) 
        addSupraUnit(); 
        
    if (col.entanglementType == toIterate || col.entanglementType == iterated) 
        addSupraUnitNeighbor(supraUnitIndexes[col.x][col.z]);
    
    if (col.entanglementType == noEntanglement) 
        supraUnitIndexes[col.x][col.z] = indexLastSupraUnit;

    
    dummy = false; // do not forget
}

void BareStateCache::pop(const Column& col)
{
    if (dummy == true) {
        dummy = false; // do not forget  
        return ; 
    }

    if (col.entanglementType == noEntanglement) 
        supraUnitIndexes[col.x][col.z] = -1; 

    if (col.entanglementType == toIterate || col.entanglementType == iterated) 
        removeSupraUnitNeighbor(supraUnitIndexes[col.x][col.z]);

    if (col.rank == 0) 
        removeLastSupraUnit(); 
  
    if (col.affectedBy() != 0) {
        for (unsigned int y = 0; y < ROWS; y++) 
            pushOrPopTritAtAAndB(false, TritAtAAndB(col, y));            
    } else {
        pushOrPopTritAtAAndB(false, TritAtAAndB(col, col.indexValue));
    }
    
    thetaEffect.addTritValue(2 * col.affectedBy(), col.x, col.z);
    parityPlane.addTritValue(2 * col.parity(), col.x, col.z);
    indexLastUnit --; 
}

void BareStateCache::iterateEntanglement(vector<Column> &unitList, Column& current)
{
    vector <bool> supraUnitFound(indexLastSupraUnit + 1 , false);
    vector<unsigned int> indexes = getUnitIndexesOfAComponent(supraUnitIndexes[current.x][current.z], supraUnitFound);

    for (vector<unsigned int>::const_iterator ind = indexes.begin(); ind != indexes.end(); ind++) {
          
        Column& col = unitList[*ind];
        col.changeScalar();
        if (col.entanglementType == noEntanglement) {
            // Donnot change the cache twice. 
            parityPlane.multiplyTritBy2(col.x, col.z); 
            thetaEffect.multiplyTritBy2(col.x, col.z);
            for (unsigned int y = 0; y < ROWS; y++) {
                TritPositionAtAAndB trit(col.x, y, col.z);
                stateA.multiplyTritBy2(trit.pA); 
                stateB.multiplyTritBy2(trit.pB); 
            }
        }
    }
    current.entanglementType = iterated; 
    current.indexValue = 0;
}

bool BareStateCache::lastSupraUnitInTheComponentOf(unsigned int index) const
{
    int supraUnit;
    list<int>::const_iterator i;
    map<int, bool> visited;
    queue<int> fifo;

    fifo.push(indexLastSupraUnit);
    visited[indexLastSupraUnit] = true;
    while (!fifo.empty()) {
        supraUnit = fifo.front(); 
        fifo.pop();
        for (i = neighborsSupraUnit[supraUnit].begin();
             i != neighborsSupraUnit[supraUnit].end(); ++i) {
            if (*i == index)
                return true;
            if (!visited[*i]) {
                fifo.push(*i); 
                visited[*i] = true; 
            }
        }
    }
    return false; 
}

vector<vector<unsigned int>> BareStateCache::getUnitIndexesOfAllTheComponents() const
{  
    vector<vector<unsigned int>> components;
    queue<int> fifo;
    vector<bool> supraUnitFound(neighborsSupraUnit.size(), false); 
    for (unsigned int indexSupraUnit = 0; indexSupraUnit != indexLastSupraUnit + 1; indexSupraUnit++) {
        if (supraUnitFound[indexSupraUnit] == false) 
            components.push_back(getUnitIndexesOfAComponent(indexSupraUnit, supraUnitFound));
    }
    return components;
}

vector<unsigned int> BareStateCache::getUnitIndexesOfAComponent(
                                unsigned int supraUnitIndexOfTheComponent, 
                                vector<bool>& supraUnitFound) const
{
    vector<unsigned int> unitIndexes; 
    queue<unsigned int> fifo;
    int supraUnit; 
    list<int>::const_iterator i;

    fifo.push(supraUnitIndexOfTheComponent);
    supraUnitFound[supraUnitIndexOfTheComponent] = true; 
   
    while (!fifo.empty()) {
        supraUnit = fifo.front(); 
        unsigned int start = startSupraUnit[supraUnit]; 
        unsigned int end   = (supraUnit == indexLastSupraUnit) ? (indexLastUnit + 1): startSupraUnit[supraUnit + 1]; 
        for (unsigned int i = start; i != end; i++) { 
            unitIndexes.push_back(i); 
        }
        fifo.pop();
        for (i = neighborsSupraUnit[supraUnit].begin();
             i != neighborsSupraUnit[supraUnit].end(); ++i) {
            if (supraUnitFound[*i] == false) {
                fifo.push(*i); 
                supraUnitFound[*i] = true; 
            }
        } 
    }
    return unitIndexes; 
}

void BareStateCache::pushOrPopTritAtAAndB(bool push, const TritAtAAndB& trit)
{
    nrActiveTrytesA -= stateA.isTheTritInAnActiveTryte(trit.pA); 
    stateA.addTritValue( (push ? 1 : 2) * trit.vA, trit.pA);
    nrActiveTrytesA += stateA.isTheTritInAnActiveTryte(trit.pA); 

    nrActiveTrytesB -= stateB.isTheTritInAnActiveTryte(trit.pB); 
    stateB.addTritValue( (push ? 1 : 2) * trit.vB, trit.pB);
    nrActiveTrytesB += stateB.isTheTritInAnActiveTryte(trit.pB);

    if (trit.stable) {
        nrStableTrytesA -= stableTritsA.isTheTritInAnActiveTryte(trit.pA); 
        stableTritsA.setTritValue( push ? (stateA.isTritActive(trit.pA)) : 0, trit.pA);
        nrStableTrytesA += stableTritsA.isTheTritInAnActiveTryte(trit.pA); 

        nrStableTrytesB -= stableTritsB.isTheTritInAnActiveTryte(trit.pB); 
        stableTritsB.setTritValue( push ? (stateB.isTritActive(trit.pB)) : 0, trit.pB);
        nrStableTrytesB += stableTritsB.isTheTritInAnActiveTryte(trit.pB);
    } else {
        if (push)
            possibleUnstableTrits.push_back(trit); 
        else 
            possibleUnstableTrits.pop_back(); 
    } 
}

void BareStateCache::addSupraUnit()
{
    neighborsSupraUnit.push_back(list<int>());
    indexLastSupraUnit++;
    startSupraUnit.push_back(indexLastUnit);
}

void BareStateCache::removeLastSupraUnit()
{
    neighborsSupraUnit.pop_back(); 
    indexLastSupraUnit--; 
    startSupraUnit.pop_back();
}

void BareStateCache::addSupraUnitNeighbor(unsigned int index)
{
    neighborsSupraUnit[index].push_back(indexLastSupraUnit); 
    neighborsSupraUnit[indexLastSupraUnit].push_back(index);
}

void BareStateCache::removeSupraUnitNeighbor(unsigned int index)
{
    neighborsSupraUnit[indexLastSupraUnit].pop_back(); 
    neighborsSupraUnit[index].pop_back(); 
} 


unsigned int TwoRoundTrailCoreCostBoundFunction::
           getCost(const vector<Column>& unitList, const BareStateCache& cache) const
{
    int contributionNewStable = 0; 
    unsigned int contributionUnstable = 0; 
    ActiveState possibleActiveTritsA = cache.stableTritsA; 
    ActiveState possibleActiveTritsB = cache.stableTritsB;
    vector<TritAtAAndB> unstable; 

    for (vector<TritAtAAndB>::const_iterator trit = cache.possibleUnstableTrits.begin(); 
            trit != cache.possibleUnstableTrits.end(); trit++) {
        if ( isTritStillUnstable(*trit, cache, unitList) == false ) {
            // The trit is now stable
            contributionNewStable -= 2 * alpha * possibleActiveTritsA.isTheTritInAnActiveTryte(trit->pA);
            if (cache.stateA.isTritActive(trit->pA))
                possibleActiveTritsA.activateTrit(trit->pA);
            contributionNewStable += 2 * alpha * possibleActiveTritsA.isTheTritInAnActiveTryte(trit->pA);

            contributionNewStable -= 2 * beta * possibleActiveTritsB.isTheTritInAnActiveTryte(trit->pB);
            if (cache.stateB.isTritActive(trit->pB))
                possibleActiveTritsB.activateTrit(trit->pB);
            contributionNewStable += 2 * beta * possibleActiveTritsB.isTheTritInAnActiveTryte(trit->pB);
            
        } else {
            // The trit is still unstable
            unstable.push_back(*trit); 
        }
    }

    // Get the contribution of the unstable trits
    for (vector<TritAtAAndB>::const_iterator trit = unstable.begin(); trit != unstable.end(); trit++) 
        contributionUnstable += getContributionUnstableTrit(*trit, possibleActiveTritsA, possibleActiveTritsB);
      
    return contributionNewStable + contributionUnstable
           + 2 * alpha * cache.nrStableTrytesA + 2 * beta * cache.nrStableTrytesB; 
}
       
bool TwoRoundTrailCoreCostBoundFunction::
           isTritStillUnstable(const TritAtAAndB& trit,
                               const BareStateCache& cache,
                               const vector<Column>& unitList) const
{
    if (cache.parityPlane.isTritActive(trit.pB.x, trit.pB.z) 
         && cache.thetaEffect.isTritActive(trit.pB.x, trit.pB.z))
            return false; 
    
    if ( trit.vA != trit. vB) {
        // The trit belongs to an affected unit. 
        if (trit.dSupraUnit < (unitList.back().dSupraUnit - 1))
            return false; 
    }
    return true; 
}

unsigned int TwoRoundTrailCoreCostBoundFunction::
           getContributionUnstableTrit(
                                 const TritAtAAndB& trit, 
                                 ActiveState& possibleActiveTritsA, 
                                 ActiveState& possibleActiveTritsB) const
{
    // Use the three activity invariants described in the paragraph 
    // "Lower bounding the cost" of Section 5.4.

    unsigned int cA, cB, contribution_1; 
    unsigned int contribution_2, contribution_3;

    // compute contribution 1
    if (possibleActiveTritsA.isTheTritInAnActiveTryte(trit.pA))
        cA = 0; 
    else
        cA = 2 * alpha;
    if (possibleActiveTritsB.isTheTritInAnActiveTryte(trit.pB))
        cB = 0; 
    else 
        cB = 2 * beta; 
    contribution_1 = min(cA, cB); 

    // contribution of unstable affected unit
    if (trit.vA != trit.vB) { 
        if (contribution_1 > 0) {
            possibleActiveTritsA.activateTrit(trit.pA); 
            possibleActiveTritsB.activateTrit(trit.pB);
            return contribution_1; 
        } else 
            return 0; 
    }

    // contribution of unstable unaffected unit
    if ( (cA == 0) && (cB == 0) )
        return 0; 
        
    // compute contribution_2 and contribution_3
    contribution_2 = 2 * alpha; 
    contribution_3 = 2 * beta; 
    for (unsigned int y = 0; y < ROWS; y++) {
        TritPositionAtAAndB pos(trit.pB.x, y, trit.pB.z);
        if (possibleActiveTritsA.isTheTritInAnActiveTryte(pos.pA))
            contribution_2 = 0; 
        if (possibleActiveTritsB.isTheTritInAnActiveTryte(pos.pB))
            contribution_3 = 0;  
    }

    if (contribution_2 > contribution_3) {
        for (unsigned int y = 0; y < ROWS; y++) {
            TritPositionAtAAndB pos(trit.pB.x, y, trit.pB.z);
            possibleActiveTritsA.activateTrit(pos.pA); 
        }
        return contribution_2; 
    }
    if (contribution_3 > contribution_2) {
        for (unsigned int y = 0; y < ROWS; y++) {
            TritPositionAtAAndB pos(trit.pB.x, y, trit.pB.z);
            possibleActiveTritsB.activateTrit(pos.pB); 
        }
        return contribution_3; 
    }
    if (contribution_1 != 0) { // contribution_3 == contribution_2 > 0
        possibleActiveTritsA.activateTrit(trit.pA); 
        possibleActiveTritsB.activateTrit(trit.pB);
        return contribution_1; 
    }
    if (contribution_2 > 0) { // and contribution_3 == 0
        for (unsigned int y = 0; y < ROWS; y++) {
            TritPositionAtAAndB pos(trit.pB.x, y, trit.pB.z);
            possibleActiveTritsA.activateTrit(pos.pA); 
        }
        return contribution_2; 
    }
    if (contribution_3 > 0) { // and contribution_2 == 0
        for (unsigned int y = 0; y < ROWS; y++) {
            TritPositionAtAAndB pos(trit.pB.x, y, trit.pB.z);
            possibleActiveTritsB.activateTrit(pos.pB); 
        }
        return contribution_3; 
    }
    return 0; 
}

void BareState::set(const vector<Column>& unitList, 
                    const BareStateCache& cache, 
                    const TwoRoundTrailCoreCostBoundFunction& costF, 
                    unsigned int aMaxCost)
{
    unsigned int cost =  costF.alpha * 2 * cache.nrActiveTrytesA
                       + costF.beta * 2 * cache.nrActiveTrytesB;
    if (cost > aMaxCost || unitList.back().ending == false) 
        valid = false; 
    else {
        valid = true;
        stateA = cache.stateA; 
        stateB = cache.stateB; 
        wA = 2 * cache.nrActiveTrytesA; 
        wB = 2 * cache.nrActiveTrytesB;

        // initialize firstActiveTritsAllowed and outKernelComponentsColumns
        firstActiveTritsAllowed = vector<UINT8>(COLUMNS * SLICES, 0x3);
        outKernelComponentsColumns.clear();

        vector<vector<unsigned int>> componentsIndexes = cache.getUnitIndexesOfAllTheComponents();
        nrComponents = componentsIndexes.size(); 
        
        for (unsigned int i = 0; i < componentsIndexes.size(); i++) {
            for (unsigned int j = 0; j < componentsIndexes[i].size(); j++) {

                unsigned int ind = componentsIndexes[i][j];
                const Column& col = unitList[ind];
                if (col.entanglementType == noEntanglement) {
                    unsigned int parity = cache.parityPlane.getTrit(col.x, col.z); 
                    unsigned int affectedBy = cache.thetaEffect.getTrit(col.x, col.z);
                    TroikaColumns troikaCol(col.x, col.z, (bool)affectedBy, parity, i);
                    if (affectedBy == 0) {
                        if (col.indexValue == 0)
                            firstActiveTritsAllowed[col.x + 9 * col.z] = 0x2;
                        if (col.indexValue == 1)
                            firstActiveTritsAllowed[col.x + 9 * col.z] = 0x4; 
                        if (col.indexValue == 2)
                            firstActiveTritsAllowed[col.x + 9 * col.z] = 0x0; 
                    } else {
                        firstActiveTritsAllowed[col.x + 9 * col.z] = 0x0; 
                        troikaCol.addValues(cache.stateB.getColumn(col.x, col.z));
                    }
                    outKernelComponentsColumns.push_back(troikaCol);
                }
            }
        }
    }
}
