/* forwardInKernelExtension.cpp */ 
#include "forwardInKernelExtension.h"
#include "state.h"

         
ForwardInKernelExtensionPreparation::ForwardInKernelExtensionPreparation
                                         (long double aMaxWeightExtension,  
                                          const ActiveState& activeB):
maxWeightExtension(aMaxWeightExtension)
{
    for (int z = 0; z < SLICES; z++) {
        for (int xTryte = 0; xTryte < 3; xTryte++) {
            for (int  y = 0; y < ROWS; y++) {
                if (activeB.isTryteActive(xTryte, y, z))
                    posActiveTrytesC.push_back(TrytePosition(xTryte, y, z)); 
            }
        }
    }
    extensionPossible = true; 
    initializePossibleActiveTritsAtCandD();
    if (extensionPossible) {
        initializeMandatoryActiveTritsAtCandD(activeB);
        int minWeightExtension = 2 * (posActiveTrytesC.size() + mandatoryActiveTritsAtC.getNrActiveTrytes());
        if ( minWeightExtension > maxWeightExtension) 
            extensionPossible = false; 
    }
}

bool ForwardInKernelExtensionPreparation::couldBeExtended() const
{
    return extensionPossible;
}

void ForwardInKernelExtensionPreparation::initializePossibleActiveTritsAtCandD()
{
    for (vector<TrytePosition>::const_iterator it = posActiveTrytesC.begin(); 
         it != posActiveTrytesC.end(); ++it) {

        for (int xOffset = 0; xOffset < 3; xOffset++) {
            TritPosition t(3 * it->x + xOffset, it->y, it->z);
            possibleActiveTritsAtC.activateTrit(t); 
            t.SRSL(); 
            possibleActiveTritsAtD.activateTrit(t); 
        }
    }

    // Filter : there cannot be only one active trit in a column at D
    for (int z = 0; z < SLICES; z++) {
        for (int x = 0; x < COLUMNS; x++) {
            UINT8 isYiActive = possibleActiveTritsAtD.getIsYiActive(x, z);
            for (int y = 0; y < ROWS; y++) {
                if (isYiActive == (1 << y)) {
                    // there is only one trit in that column
                    TritPosition t(x, y, z); 
                    possibleActiveTritsAtD.deactivateTrit(t); 
                    t.invSRSL(); 
                    possibleActiveTritsAtC.deactivateTrit(t); 
                    if (!possibleActiveTritsAtC.isTryteActive(t)) {
                        extensionPossible = false; 
                        return; 
                    }
                }
            }
        }
    }
}

void ForwardInKernelExtensionPreparation::initializeMandatoryActiveTritsAtCandD(const ActiveState& activeB)
{
    int activeTryteAtB; 
    int activeTryteAtC; 
    TritPosition t;

    for (vector<TrytePosition>::const_iterator it = posActiveTrytesC.begin(); 
         it != posActiveTrytesC.end(); ++it) {
        // with the properties of the S-box
        // activeTryteAtB == 0x4 : | 0 | 0 | X | --> | . | . | X |
        // activeTryteAtB == 0x1 : | X | 0 | 0 | --> | . | X | . |                
        // if there is only one active trit in an active tryte at C, 
        // this trit is mandatory
        activeTryteAtB = activeB.getActiveTryte(it->x, it->y, it->z); 
        if (activeTryteAtB == 0x1 || activeTryteAtB == 0x4) {
            if (activeTryteAtB == 0x4)
                t.set(3 * it->x + 2, it->y, it->z); 
            else 
                t.set(3 * it->x + 1, it->y, it->z); 
            if (!possibleActiveTritsAtC.isTritActive(t)) {
                extensionPossible = false; 
                return; 
            }
            mandatoryActiveTritsAtC.activateTrit(t); 
            t.SRSL(); 
            mandatoryActiveTritsAtD.activateTrit(t);   
            addIfNeededATwoTritsMandatoryColumnAtD(t.x, t.z); 
        }
        // if there is only one active trit in an active tryte at C, 
        // this trit is mandatory
        activeTryteAtC = possibleActiveTritsAtC.getActiveTryte(it->x, it->y, it->z);
        for (int tritIndex = 0; tritIndex < 3; tritIndex++) {
            if (activeTryteAtC == ( 1 << tritIndex )) {
                t.set(3 * it->x + tritIndex, it->y, it->z);
                mandatoryActiveTritsAtC.activateTrit(t); 
                t.SRSL(); 
                mandatoryActiveTritsAtD.activateTrit(t); 
                addIfNeededATwoTritsMandatoryColumnAtD(t.x, t.z); 
            }
        }
    }
    
}

void ForwardInKernelExtensionPreparation::addIfNeededATwoTritsMandatoryColumnAtD(int x, int z)
{
    if (possibleActiveTritsAtD.getIsYiActive(x, z) != 0x7) {
        // this column has 2 possible active trits. When this function is 
        // called, we know that one of this tris is mandatory. So both 
        // trits are mandatory. 
        for (int y = 0; y < ROWS; y++) {
            TritPosition t(x, y, z); 
            if (possibleActiveTritsAtD.isTritActive(t)) {
                mandatoryActiveTritsAtD.activateTrit(t); 
                t.invSRSL(); 
                mandatoryActiveTritsAtC.activateTrit(t); 
            }
        }
    }
}

ostream & operator << (ostream& fout, const ForwardInKernelExtensionPreparation& prep)
{
    fout << "possible: " << prep.couldBeExtended() << endl; 
    fout << "possibleActiveTritsAtC: " << endl << prep.possibleActiveTritsAtC << endl; 
    fout << "possibleActiveTritsAtD: " << endl << prep.possibleActiveTritsAtD << endl;
    fout << "mandatoryActiveTritsAtC: " << endl << prep.mandatoryActiveTritsAtC << endl; 
    fout << "mandatoryActiveTritsAtD: " << endl << prep.mandatoryActiveTritsAtD << endl;
    return fout; 
}    

ForwardInKernelExtensionIterator::ForwardInKernelExtensionIterator(
                const ForwardInKernelExtensionPreparation& prep,
                const TroikaState& stateB)
    :maxWeightExtension(prep.maxWeightExtension), 
     possibleStatesD(prep.possibleActiveTritsAtD, prep.mandatoryActiveTritsAtD), 
     extension(stateB, prep.posActiveTrytesC)
{
    end = false;
    for (; !possibleStatesD.isEnd(); ++possibleStatesD) { 
        extension.setStateCAndDFromStateD(*possibleStatesD);
        if (extension.isValidAndBelowWeight(maxWeightExtension))
            return; 
    }
    end = true; 
}

bool ForwardInKernelExtensionIterator::isEnd() const
{
    return end; 
}

void ForwardInKernelExtensionIterator::operator++()
{
    while (!possibleStatesD.isEnd()) {
        ++possibleStatesD;
        extension.setStateCAndDFromStateD(*possibleStatesD);
        if (extension.isValidAndBelowWeight(maxWeightExtension))
            return; 
    }
    end = true; 
}

const ForwardExtension& ForwardInKernelExtensionIterator::operator*() const
{
    return extension; 
}
