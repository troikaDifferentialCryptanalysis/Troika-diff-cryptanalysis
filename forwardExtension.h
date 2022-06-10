/* forwardExtension.h */
#ifndef FORWARD_EXTENSION_H 
#define FORWARD_EXTENSION_H

/* 
  This code uses functions and ideas from KeccakTools
  (https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
  and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
  We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

  This file contains the classes used to extend a trail core in the forward
  direction(see Appendix C.2). For a forward extension (C, D), used to extend
  a trail core (..., B), the state  ShiftLanesoShiftRows(C), abbreviated 
  SRSLC, is chosen trit by trit, respecting some constraints. All the trits of
  a slice of SRSLC must be chosen before choosing the trits of another slice. 
  The weight of the extension is computed while choosing the slices. 
*/


#include <map>
#include <stack>
#include "types.h"
#include "state.h"
#include "sbox.h"
#include "extensionsIterator.h"

enum ConstraintForTritValue {noConstraint, mustBe1, mustBe2, cannotBe0};

typedef GenericExtensionIterator< 
        class ForwardExtensionPreparation, 
        class TritInfo,
        class ForwardExtensionCache,
        class CostFunctionForwardExtension,
        class ForwardExtension> 
        ForwardExtensionIterator;

class ForwardExtensionCache;

/** This class is used to store the characteristics of a trit of SRSLC 
  * and to choose its value.
  */ 
class TritInfo {
    public:  
        TritPosition posAtSRSLC;
        TritPosition posAtC; 
        /** true if the trit is the last possible active trit of its slice,
          * false otherwise.
          */
        bool mustCalculateTheCostOfTheSlice;
        /** true if  : 
          * 1) the trit is the last possible active trit of its slice of 
          *    coordinate z
          * 2) the slice of coordinate (z - 1) of SRSL(C) is passive
          * false otherwise.
          */
        bool mustCalculateTheCostOfTheNextSlice;
        unsigned int value; 
    public:
        /** The constructor. It initializes the attributes posAtSRSLC,
          * posAtC, and value. The attributes mustCalculateTheCostOfTheSlice
          * and mustCalculateTheCostOfTheNextSlice must be initialized later.
          * @param x The x coordinate of the trit at SRSL(C). 
          * @param y The y coordinate of the trit at SRSL(C). 
          * @param z The z coordinate of the trit at SRSL(C). 
          */ 
        TritInfo(unsigned int x, unsigned int y, unsigned int z); 
        void setFirstValue(const ForwardExtensionCache& cache);
        bool setNextValue(const ForwardExtensionCache& cache);
        friend ostream & operator << (ostream &fout, const TritInfo &aTritInfo);
    private: 
        /** It checks wheter the current value for the trit is valid or not.
          * That trit must respect the constraints stored 
          * in constraintAtC. Furthermore, if the value chosen for 
          * the trit is 0, the method checks if that does not deactivate a tryte
          * of C which has to be active. 
          * @param constraintAtC   Map between an active trit position
          *                        of C and the constraint for the 
          *                        value of the trit. 
          * @param nrNonActiveTritsAtC  Map between an active tryte position 
          *                             of C and the number of passive
          *                             trits that the tryte contains. 
          * @return true if the trit can take the value @value.
          */
        bool isAValidTritValue(
                const map<TritPosition, ConstraintForTritValue>& constraintAtC,
                const map<TrytePosition, unsigned int> &nrNonActiveTritsAtC) const;
};

/** Given the state B, the class initializes the attributes needed for 
  * the extension. 
  */ 
class ForwardExtensionPreparation 
{
    public:
        /** Maximum weight for w(B--ST-->C) + wMinDir(D). */ 
        long double maxWeightExtension; 
        /** Vector that stores, for all the (3 * nr of active trytes at B) 
          * possible active trits of SRSL(C), information needed about that trit
          * in order to choose its value.
          */ 
        vector<TritInfo> tritsInfoAtSRSLC; 
        /** Vector which stores all the active trytes'positions of the state B. */
        vector <TrytePosition> posForSTCompatibility;
        /** Map between an active trit position of C and the constraint for the 
          * value of the trit. 
          */
        map<TritPosition, ConstraintForTritValue> constraintAtC;
    private: 
        ActiveState possibleActiveTritsAtSRSLC;
    public: 
        ForwardExtensionPreparation(const TroikaState& stateB, 
                                    long double maxWeightExtension); 
        vector<TritInfo>& getPartsList(){ return tritsInfoAtSRSLC; }
        bool couldBeExtended() {return true; } 
        friend ostream & operator << (ostream& fout, const ForwardExtensionPreparation& prep);
    private: 
        void initPosForSTCompatibility(const TroikaState& stateB);
        void initPossibleActiveTritsAtSRSLCAndConstraintForActiveTritsAtC(const TroikaState& stateB);
        void addInfoOfTheTritsOfTheSlice(unsigned int z);
}; 

/** Auxiliary class for the class ForwardExtensionIterator.
  * It stores information to compute the cost of the extension that is 
  * being chosen.
  */ 
class ForwardExtensionCache 
{
    public: 
        /** Map between an active trit position of C and the constraint for the 
          * value of the trit. 
          */
        map<TritPosition, ConstraintForTritValue>& constraintAtC;
        /** Map between an active tryte position of C and the number of passive
          * trits that the tryte contains. 
          */
        map<TrytePosition, unsigned int> nrNonActiveTritsAtC;
        /** Trit activity pattern of the state D. */ 
        ActiveState activeD; 
        /** Each time an active slice of SRSL(C) is chosen, the accumulated
          * number of trytes of D is pushed to the stack. 
          */
        stack<unsigned int> stackNrActiveTrytesD;  
        /** For  0 ≤ x < 9,  0 ≤ y < 3,  0 ≤ z < 27, tritsAtSRSLC[x][y][z] 
          * contains the value of the trit of coordinates x, y and z of the
          * state SRSL(C). This array is updated each time a new trit value of
          * SRSL(C) is chosen.
          */ 
        unsigned int tritsAtSRSLC[COLUMNS][ROWS][SLICES];
        /** For  0 ≤ x < 9,  0 ≤ z < 27, columnParityOfSRSLC[x][z] contains 
          * the parity of the column of SRSL(C) of coordinates x and z. This 
          * array is updated each time a new slice of SRSL(C) is chosen.
          */
        unsigned int columnParityOfSRSLC[COLUMNS][SLICES];
        /** For  0 ≤ x < 9,  0 ≤ y < 3,  0 ≤ z < 27, tritsAtD[x][y][z] contains
          * the value of the trit of coordinates x, y and z of state D.
          * This array is updated each time a new slice of SRSL(C) is chosen.
          */ 
        unsigned int tritsAtD[COLUMNS][ROWS][SLICES];
    public: 
        ForwardExtensionCache(ForwardExtensionPreparation& prep, const TroikaState& stateB);
        friend ostream & operator << (ostream& fout, const ForwardExtensionCache& cache);
        /** It updates the array @a tritsAtSRSLC with the new trit value.
          * If the trit value is 0, the map nrNonActiveTritsAtC is also updated. 
          * If the trit is the last trit of it slice at SRSL(C) that must be 
          * chosen, it also updates @stackNrActiveTrytesAtD, using the arrays 
          * tritsAtSRSLC, columnParityOfSRSLC and tritsAtD. 
          * @param trit The trit of SRSL(C) which has just been chosen. 
          */
        void push(const TritInfo &trit); 
        /** It removes the contribution of @a trit to the cost : 
          * it updates stackNrActiveTrytesAtD. If the trit value is 0, the map 
          * nrNonActiveTritsAtC is also updated. 
          * @param trit The trit of SRSL(C) which is going to be removed. 
          */
        void pop(const TritInfo &trit);
    private: 
        void setTheColumnParityOfTheSliceAtSRSLC(unsigned int z);
        void applyToTheSliceAddColumnParity(unsigned int z);
        void addActiveTritsOfTheSliceAtD(unsigned int z) const;
        void removeActiveTritsOfTheSliceAtD(unsigned int z) const;
        unsigned int getNrActiveTrytesOfTheSliceAtD(unsigned int z) const;
};

/** Class used by the ForwardExtensionIterator to test if the extension 
  * that is being formed or is already formed has a too high weight.
  */ 
class CostFunctionForwardExtension
{
    public:
        /** Maximum weight for w(B-->C) + wMinDir(D). */ 
        long double maxWeightExtension; 
        unsigned int nrActiveTrytesC;
    public: 
        CostFunctionForwardExtension(ForwardExtensionPreparation &prep, const TroikaState& stateB); 
        bool tooHighCost(const ForwardExtensionCache& cache, int indCurPart) const;

}; 

/** The output representation for a forward extension (C, D) of a trail core 
  * ( ..., B ).
  */
class ForwardExtension {
    /* A k-rounds trail core  --Lambda--> B is extended into a 
     * k+1-rounds trail cores --Lambda--> B --ST--> C --L--> D.
     */
    public: 
        const TroikaState stateB;
        TroikaState stateC; 
        TroikaState stateD; 
        /** weight w(B--ST-->C). */ 
        Weight wBC;
        /** Minimum direct weight of state C. */ 
        Weight wMinDirD;
        Sbox sbox;
        /** Vector that stores the position of the active trytes of the state B. */ 
        vector<TrytePosition> posForSTCompatibility;
        /** True if stateB is compatible through Subtrytes with stateC. */ 
        bool valid; 
    public: 
        ForwardExtension(const ForwardExtensionPreparation& prep, const TroikaState &stateB);
        ForwardExtension(const TroikaState& aStateB, const vector<TrytePosition> posForSTCompatibility);
        /** It returns true if the extension is valid and has a weight below
          * @a maxWeightExtension.
          */ 
        bool isValidAndBelowWeight(double long maxExtensionWeight) const;
        void setStateCAndDFromStateD(const TroikaState& aStateD); 
        /** This method is called when all the parts of the extension are chosen.
          * It updates the attributes of the ForwardExtension using information 
          * stored in @tritsInfoAtSRSLC and @cache. 
          */ 
        void set(const vector<TritInfo>& tritsInfoAtSRSLC,
                 const ForwardExtensionCache& cache, 
                 const CostFunctionForwardExtension& costF);
    private: 
        Weight getWeight() const{return wBC + wMinDirD;}
        Weight getNrActiveTrytesC() const;        
};

template<>                     
bool ForwardExtensionIterator::next();

#endif
