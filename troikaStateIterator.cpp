/** troikaStateIterator.cpp */ 
#include "troikaStateIterator.h"
#include "mixedStateIterator.h"

TroikaColumns::TroikaColumns(unsigned int aX, unsigned int aZ, bool affected, 
                             unsigned int aParity, int componentIndex)
{
    x = aX; 
    z = aZ;
    isAffected = affected;
    parity = aParity; 
    indexOutKernelComponent = componentIndex;
    indexValue = 0; 
}

bool TroikaColumns::next()
{
    if (indexValue == columnValues.size() - 1)
        return false; 
    indexValue++; 
    return true;

}

unsigned int TroikaColumns::getTrit(unsigned int y) const
{
    return columnValues[indexValue][y]; 
}

void TroikaColumns::setValues(const vector<array<int, 3>>& values)
{
    columnValues = values; 
}

void TroikaColumns::addValues(const array<int, 3>& value)
{
    columnValues.push_back(value); 
}

ostream & operator << (ostream &fout, const TroikaColumns &aColumn)
{
    fout << (ColumnPosition)aColumn; 
    fout << " index : " << aColumn.indexValue; 
    fout << " component : " << aColumn.indexOutKernelComponent;
    fout << " nr columns: " << aColumn.columnValues.size();
    fout << " parity: " << aColumn.parity; 
    fout << " affected: " << aColumn.isAffected; 
    fout << " component: " << aColumn.indexOutKernelComponent; 
    return fout; 
}

static vector <array<int, 3>> **getPossibleColumnValues() {
    unsigned int sumsTo; 
    unsigned int isYiActive; 

    vector <array<int, 3>> **possibleColumnValues = 
                                             new vector <array<int, 3>> *[3];
    for(int i = 0; i < 3; i++) {
        possibleColumnValues[i] = new vector<array<int, 3>>[8];
    }

    for (int trit_0 = 0; trit_0 < 3; trit_0++) {
        for (int trit_1 = 0; trit_1 < 3; trit_1++) {
            for (int trit_2 = 0; trit_2 < 3; trit_2++) {
                sumsTo = (trit_0 + trit_1 + trit_2) % 3;
                array<int, 3> aColumnValue = {trit_0, trit_1, trit_2}; 
                
                isYiActive = 0; 
                for (unsigned int y = 0; y < 3; y++)
                    if (aColumnValue[y])
                        isYiActive += (1 << y); 
                possibleColumnValues[sumsTo][isYiActive].push_back(aColumnValue); 
            }
        }
    }
    return possibleColumnValues; 
}

vector <array<int, 3>> **TroikaStateIterator::possibleColumnValues
                                                  = getPossibleColumnValues();

const array<int, 3> TroikaStateIterator::inKernelColumn[9] = { {0, 0, 0},
                                                               {1, 2, 0}, 
                                                               {2, 1, 0}, 
                                                               {1, 0, 2}, 
                                                               {2, 0, 1}, 
                                                               {0, 1, 2}, 
                                                               {0, 2, 1},
                                                               {1, 1, 1}, 
                                                               {2, 2, 2}};

TroikaStateIterator::TroikaStateIterator(const ActiveState& activeInKernel)
: nrOutKernelComponents(0), multiplyOutKernelComponent(0)
{
    unsigned int isYiActive; 
    for (unsigned int z = 0; z < SLICES; z++) {
        for (unsigned int x = 0; x < COLUMNS; x++) {
            isYiActive = activeInKernel.getIsYiActive(x, z); 
            if (isYiActive) {
                TroikaColumns columns(x, z);
                columns.setValues(possibleColumnValues[0][isYiActive]); 
                columnsList.push_back(columns);
            }
        }
    }
    if (columnsList.size() == 0) {
        end = true; 
    } else {
        end = false; 
        first(); 
    } 

}

TroikaStateIterator::TroikaStateIterator(const ActiveState& possibleActiveTrits, 
                                         const ActiveState& mandatoryActiveTrits)
:nrOutKernelComponents(0), multiplyOutKernelComponent(0)
{
    bool isAValidColumn; 
    for (unsigned int z = 0; z < SLICES; z++) {
        for (unsigned int x = 0; x < COLUMNS; x++) {
            if (possibleActiveTrits.isColumnActive(x, z)) {
                TroikaColumns columns(x, z);

                // search in the array inKernelColumn the valid column values
                for (unsigned int i = 0; i < NR_IN_KERNEL_COLUMN; i++) {
                    isAValidColumn = true;
                    // apply to the 3 trits of inKernelColumn[i] a filter to
                    // know if inKernelColumn[i] is a valid column value
                    for (unsigned int y = 0; y < ROWS; y++) { 
                        if (inKernelColumn[i][y]) {
                            if (!possibleActiveTrits.isTritActive(x, y, z))
                                isAValidColumn = false; 
                        } else {
                            if (mandatoryActiveTrits.isTritActive(x, y, z))
                                isAValidColumn = false; 
                        }
                    }
                    if (isAValidColumn) 
                        columns.addValues(inKernelColumn[i]);
                }
                columnsList.push_back(columns); 
            }
        }
    }
    if (columnsList.size() == 0) {
        end = true; 
    } else {
        end = false; 
        first(); 
    } 
}

TroikaStateIterator::TroikaStateIterator(const ActiveTrailCoreCache& cache)
: nrOutKernelComponents(0), multiplyOutKernelComponent(0)
{
    vector<TroikaColumns>::const_iterator it; 
    unsigned int isYiActive; 
    for (it = cache.inKernelColumns.begin(); it != cache.inKernelColumns.end(); it++) {
        isYiActive = cache.activeB.getIsYiActive(it->x, it->z);
        TroikaColumns col = *it; 
        col.setValues(possibleColumnValues[0][isYiActive]);  
        columnsList.push_back(col); 
    }
    if (columnsList.size() == 0) {
        end = true; 
    } else {
        end = false; 
        first(); 
    }
}

TroikaStateIterator::TroikaStateIterator(const MixedTrailCoreCache& cache)
:nrOutKernelComponents(cache.nrComponents), multiplyOutKernelComponent(0)
{
    vector<TroikaColumns>::const_iterator it;
    unsigned int isYiActive;
    // Add the columns that do not belong to a component, ie the unaffectec
    // columns of zero parity. 
    for (it = cache.inKernelColumns.begin(); it != cache.inKernelColumns.end(); it++) {
        TroikaColumns col = *it; 
        isYiActive = cache.activeB.getIsYiActive(it->x, it->z);
        col.setValues(possibleColumnValues[0][isYiActive]);  
        columnsList.push_back(col); 
    }
    // Add the columns that belong to a component.
    for (it = cache.outKernelComponentColumns.begin(); it != cache.outKernelComponentColumns.end(); it++) {
        TroikaColumns col = *it;
        isYiActive = cache.activeB.getIsYiActive(col.x, col.z);
        // If the column is an unaffected column of a run, its attribute columnValues
        // has to be initialized, otherwise the column is an affected column and 
        // its attribute columnValues is already initialized. 
        if (col.isAffected == false) // unaffected column of a run  
            col.setValues(possibleColumnValues[col.parity][isYiActive]);
        columnsList.push_back(col); 
    }
    if (columnsList.size() == 0) {
        end = true; 
    } else {
        end = false; 
        first(); 
    }
}

bool TroikaStateIterator::isEnd() const
{
    return end; 
}

void TroikaStateIterator::operator++()
{
    if (end)
        return;

    if (nrOutKernelComponents > 0) {
        if ((int)multiplyOutKernelComponent < (1 << nrOutKernelComponents) - 1) {
            multiplyOutKernelComponent ++; 
            return;
        }
        multiplyOutKernelComponent = 0; 
    }
    do {
        if (toSibling())
            break; 

        if (!toParent()) {     
            end = true; 
            return;
        } 
    } while (true);
 
    do {
        if (!toChild())
            return; 

    } while (true);
}

const TroikaState& TroikaStateIterator::operator*()
{
    int tritValue; 
    vector<TroikaColumns>::const_iterator it; 

    state.setEmptyState(); 
    for (it = columnsList.begin(); it != columnsList.end(); it++) {
        for (unsigned int y = 0; y < 3; y++) {
            if (multiplyOutKernelComponent & (1 << it->indexOutKernelComponent))
                tritValue = (2 * (*it).getTrit(y)) % 3; 
            else
                tritValue = (*it).getTrit(y);
            state.setTritValue(tritValue, it->x, y, it->z);
        } 
    }

    return state; 

}

void TroikaStateIterator::first()
{
    vector<TroikaColumns>::iterator it; 
    for (it = columnsList.begin(); it != columnsList.end(); it++)
        it->indexValue = 0; 
    indexColumn = columnsList.size() - 1; 
}

bool TroikaStateIterator::toChild()
{
   if (indexColumn == (int)columnsList.size() - 1)
        return false; 
   else {
        indexColumn++; 
        columnsList[indexColumn].indexValue = 0; 
        return true; 
    }
}

bool TroikaStateIterator::toSibling()
{
    return columnsList[indexColumn].next(); 
}

bool TroikaStateIterator::toParent()
{
    if (indexColumn == 0)
        return false; 
    else {
        indexColumn--;
        return true; 
    }
}  

ostream & operator << (ostream &fout, const TroikaStateIterator &it)
{
    fout << "end: " << it.isEnd() << endl;
    fout << "nrOutKernelComponents: " << it.nrOutKernelComponents << endl; 
    fout << "multiplyOutKernelComponent: " << it.multiplyOutKernelComponent << endl; 
    fout << "indexColumn: " << it.indexColumn << endl;
    fout << "COLUMNS" << endl; 
    for (unsigned int i = 0; i < it.columnsList.size(); i++)
        fout << it.columnsList[i] << endl; 
    return fout;
}
