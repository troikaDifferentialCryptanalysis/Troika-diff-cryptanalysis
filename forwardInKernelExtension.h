/* forwardInKernelExtension.h */
#ifndef FORWARD_IN_KERNEL_EXTENSION_ITERATOR_H 
#define FORWARD_IN_KERNEL_EXTENSION_ITERATOR_H 

/*
    This code uses functions and ideas from KeccakTools
    (https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
    and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
    We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

   The classe of this file are used to extend a trail core in the forward 
   direction and in the kernel (see Appendix C.3)
*/

#include "state.h"
#include "mixedStateIterator.h"
#include "trailCore.h"
#include "forwardExtension.h"
#include "troikaStateIterator.h"

/** Class used to see if trail core ( ..., B ) defined by the trit activity 
  * pattern of B can possibly be extended into a trail core 
  *  ... --> B --ST--> C --Lambda--> D with D in the kernel. The algorithm is 
  *  described in Appendix C.3.
  */ 
class ForwardInKernelExtensionPreparation {
    public:
        /* Maximum weight for w(B-->C) + wMinDir(D). */ 
        const long double maxWeightExtension; 
        vector <TrytePosition> posActiveTrytesC;
        ActiveState possibleActiveTritsAtC; 
        ActiveState mandatoryActiveTritsAtC; 
        ActiveState possibleActiveTritsAtD; 
        ActiveState mandatoryActiveTritsAtD;
    private: 
        /** True if the 2 tests that remove some impossible extensions
          * were successful. 
          */
        bool extensionPossible;
    public:
        /* The constructor. It initilizes all the attributes. 
         * @param aMaxWeightExtension maximum for the weight w(B-->C) + wMinDir(D)
         * @param activeB The trit activity pattern of B
         */
        ForwardInKernelExtensionPreparation(long double aMaxWeightExtension,  
                                            const ActiveState& activeB);
        bool couldBeExtended() const; 
        friend ostream & operator << (ostream &fout, const ForwardInKernelExtensionPreparation& prep);
    private: 
        /** It is the first filter. It initializes the variables 
          * possibleActiveTritsAtC and possibleActiveTritsAtD. Each active trytes
          * of B indicates 3 possible active trits for C / D. Before the filter,
          * the active states possibleActiveTritsAtC and possibleActiveTritsAtD 
          * have (3 * number active trytes of B) active trits. Then, the active
          * trits of D which are alone on their column are deactivated as well 
          * as the corresponding trits of C. We check if it doesn't deactivate a
          * tryte of C necessarily active.
          */
        void initializePossibleActiveTritsAtCandD();
        /** It is the second filter. It initializes the variables
          * mandatoryActiveTritsAtC and mandatoryActiveTritsAtD using properties
          * of the Sbox and the localisation of the active trytes of C (at the 
          * same positions as those of D).
          */
        void initializeMandatoryActiveTritsAtCandD(const ActiveState& activeB);
        /** If the column of D of coordinates @x an @z has exactly 2 possible
          * active trits, including one that is necessarily active, the second
          * is necessarily active. In that case, adds to the variable 
          * mandatoryActiveTritsAtC and mandatoryActiveTritsAtD the positions
          * of the trits necessarily active.
          */
        void addIfNeededATwoTritsMandatoryColumnAtD(int x, int z);

}; 

/** Class used to find all the valid forward in kernel extensions. */ 
class ForwardInKernelExtensionIterator {
    private: 
        /** The maximum weight for w(B--ST-->C) + wMinDir(D). */ 
        const long double maxWeightExtension; 
        /** Variable used to store a valid extension. */ 
        ForwardExtension extension;
        /** Variable used to find all the states D that are in the kernel and 
          * that respect the constaints given by the variables  
          * possibleActiveTritsAtD and mandatoryActiveTritsAtD of the 
          * PreparationForForwardInKernelExtension.
          */ 
        TroikaStateIterator possibleStatesD; 
        bool end;
    public:
        /** The constructor. @prep is a PreparationForForwardInKernelExtension
          * that has passed the test. If it is possible to extend the trail core
          * (..., B), the first extension (C, D) found is used to initialize 
          * the attribute extension. 
          * @param aStateAtB   The state 
          * @param preparation A PreparationForForwardInKernelExtension
          *                    associated to a trail core (A, B) which passes 
          *                    the test.
          */
        ForwardInKernelExtensionIterator(
                const ForwardInKernelExtensionPreparation& prep, 
                const TroikaState& stateB);
        /** It indicates if the last extension has been reached. */
        bool isEnd() const;
        /** The '++' operator moves the iterator to the next extension. */
        void operator++();
        /** It returns a constant reference to the current extension. */
        const ForwardExtension& operator*() const;
};

#endif
