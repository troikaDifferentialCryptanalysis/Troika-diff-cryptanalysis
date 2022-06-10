/* backwardExtension.h */
#ifndef BACKWARD_EXTENSION_H 
#define BACKWARD_EXTENSION_H 

/* 
  This code uses functions and ideas from KeccakTools
  (https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
  and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
  We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

  This file contains the classes used to extend a trail core in the backward
  direction(see Appendix C.2). For a backward extension (A, B), used to extend
  a trail core (C, ...), the state B is chosen active tryte by active tryte.
  All the trytes of a slice of B must be chosen before choosing the trytes of 
  another slice. The weight of the extension is computed while choosing the slices. 
*/

#include <stack>
#include "types.h"
#include "state.h"
#include "sbox.h"
#include "extensionsIterator.h"

typedef GenericExtensionIterator< 
        class BackwardExtensionPreparation, 
        class TryteInfo,
        class BackwardExtensionCache,
        class CostFunctionBackwardExtension,
        class BackwardExtension> 
        BackwardExtensionIterator;

class BackwardExtensionCache; 

/** This class is used to store all the possible values that an active tryte
  * of B can take. 
  */ 
class TryteInfo {
    public: 
        /** The position of the tryte. */
        TrytePosition position; 
        /** Vector that stores all the possible values for the tryte. */
        vector<TryteSTCompatible> possibleValues; 
        /** possibleValues[index] is the current tryte value. */
        unsigned int index; 
        /** true if the tryte is the last active tryte of its slice,
          * false otherwise.
          */
        bool mustCalculateTheCostOfTheSlice;
        /** true if : 
          * 1) the tryte is the last active tryte of its slice of 
          *    coordinate z
          * 2) the slice of B of coordinate (z - 1)  is passive
          * false otherwise.
          */
        bool mustCalculateTheCostOfTheNextSlice;
    public:

        /* The constructor. The attributes  mustCalculateTheCostOfTheSlice 
         * and mustCalculateTheCostOfTheNextSlice must be initialized later. 
         * @param aPosition The position of the active tryte.  
         * @parame tryteAtC The tryte of C that has the same position. 
         *                  The values stored in possibleValues must be 
         *                  compatible through Subtrytes with tryteAtC.
         */
        TryteInfo(const TrytePosition& aPosition, const Tryte & tryteAtC);
        /** It sets  the first tryte value.
          * @param cache Not useful here.
          */ 
        void setFirstValue(const BackwardExtensionCache& cache);
        /** It sets (if possible) the next value of the tryte. 
          * @return true if another value has been set, false otherwise.
          */
        bool setNextValue(const BackwardExtensionCache& cache);
        /** @return The current tryte value. */
        Tryte getValue() const; 
        /** @return The current weight of the transion. */
        Weight getWeight() const;
        friend ostream & operator << (ostream &fout, const TryteInfo &aTryteInfo);
};

/** Given a state C, this class is used to initialize the vector partsList 
  * of the BackwardExtensionIterator.
  */ 
class BackwardExtensionPreparation 
{
    public:
        long double maxWeightExtension;
        vector<TryteInfo> trytesInfoAtB;
    public: 
        BackwardExtensionPreparation(const TroikaState &stateC, long double maxWeightExtension); 
        /* @return A reference to the vector trytesInfoAtB. */
        vector<TryteInfo>& getPartsList();
        bool couldBeExtended() {return true; }
    private:
        void addInfoOfTheTrytesOfTheSliceAtB(unsigned int z,
                                             const TroikaState &stateAtC); 
        friend ostream & operator << (ostream &fout, const BackwardExtensionPreparation& prep);
};

/** Auxiliary class for the class BackwardExtensionIterator.
  * It stores information to compute the cost of the extension.
  */ 
class BackwardExtensionCache 
{
    public: 
        /** Each time an active tryte of B is chosen, the accumulated weight 
          * w(B --ST-->C) is updated and the new value is pushed to the stack.
          * That weight is the sum of all the weights of the trytes of B that
          * were chosen. 
          */
        stack<Weight> stackWeightBC;
        /** The trit activity pattern of the state A of the extension. */ 
        ActiveState activeA;
        /** The number of active trytes of the state activeA. */ 
        unsigned int nrActiveTrytesA; 
        /** For  0 ≤ x < 9,  0 ≤ y < 3,  0 ≤ z < 27, 
          * B[x][y][z] contains the value of the trit of stateB of coordinates
          * x, y and z.
          */ 
        unsigned int B[COLUMNS][ROWS][SLICES];
        /** For  0 ≤ x < 9,  0 ≤ y < 3,  0 ≤ z < 27, 
          * invAddColumnParityOfB[x][y][z] contains the value of the trit 
          * of coordinates x, y and z of AddColumnParity^{-1} (B).
          */ 
        unsigned int invAddColumnParityOfB[COLUMNS][ROWS][SLICES];
        /** For  0 ≤ x < 9,  0 ≤ z < 27, columnParityOfB[x][z] contains the
          * parity of the column of B coordinates x and z. 
          */
        unsigned int columnParityOfB[COLUMNS][SLICES];
    public:
        BackwardExtensionCache(const BackwardExtensionPreparation &prep, 
                               const TroikaState &stateC);
        /** It updates the array @a B with the new tryte value and the 
          * stackWeightBC with the new weight value. If the tryte is 
          * the last tryte of its slice that must be chosen, it also updates 
          * @trytesAtA, using the arrays columnParityOfB, invAddColumnParityOfB
          * and B.
          * @param tryteInfoAtB The tryte of B which has just been chosen. 
          */
        void push(const TryteInfo& tryteAtB);
        /** It removes the contribution of tryteInfoAtB to the cost : 
          * it updates stackWeightBC and trytesAtA.
          * @param tryteInfoAtB The tryte of B which is going to be removed. 
          */
        void pop(const TryteInfo& tryteAtB);
        friend ostream & operator << (ostream &fout, const BackwardExtensionCache& cache);
    private: 
        void initTryteOfStateB(const TryteInfo& tryteAtB);
        void applyToTheSliceInvAddColumnParity(unsigned int z);
        void setTheColumnParityOfTheSlice(unsigned int z); 
        void addTrytesAtAFromASliceOfInvAddColumnParityOfB(unsigned int z); 
        void removeTrytesAtAFromASliceOfInvAddColumnParityOfB(unsigned int z);
}; 

/** Class used by the BackwardExtensionIterator to test if the extension 
  * that is being formed or is already formed has a too high weight.
  */ 
class CostFunctionBackwardExtension
{
    public: 
        long double maxWeightExtension;
        unsigned int nrActiveTrytesB;

    public: 
        CostFunctionBackwardExtension(const BackwardExtensionPreparation& prep, 
                                      const TroikaState& stateC);
        bool tooHighCost(const BackwardExtensionCache& cache,
                          int indCurPart) const; 
};

/** The output representation for a backward extension. */
class BackwardExtension
{
    /* A k-rounds trail core C --L--> D ... --L--> ... is extended into a 
     * k+1 rounds trail cores A --L--> B --ST--> C --L--> D ... --L--> ...
     */
    public: 
        /** It stores the state A of the backward extension. */ 
        TroikaState stateA;
        /** It stores the state B of the backward extension. */ 
        TroikaState stateB; 
        /** It stores the minimum reverse weight of the state A. */
        Weight wMinRevA; 
        /** It stores the weight w(B --ST--> C) of the backward extension. */
        Weight wBC;
    public: 
        BackwardExtension():wMinRevA(0), wBC(0) {}
        BackwardExtension(BackwardExtensionPreparation& prep, const TroikaState& stateC) {} 
        /** This method is called when all the parts of the extension are chosen.
          * It updates the attributes of the BackwardExtension using information 
          * stored in @trytesInfoAtB and @cache. 
          */ 
        void set(const vector<TryteInfo> trytesInfoAtB,
                 const BackwardExtensionCache& cache, 
                 const CostFunctionBackwardExtension& costF);
        /** By construction, stateB is compatible through Subtrytes with 
          * the state C. It returns true if the weight of the extension is 
          * less than maxWeightExtension.
          */ 
        bool isValidAndBelowWeight(long double maxWeightExtension) const;
    private: 
        Weight getWeight() const{return wBC + wMinRevA;}
};

#endif 
