/* backwardInKernelExtension.h */
#ifndef BACWARD_IN_KERNEL_EXTENSION_H 
#define BACWARD_IN_KERNEL_EXTENSION_H 

/* 
   This code uses functions and ideas from KeccakTools
   (https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
   and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
   We thank Silvia Mella and Gilles Van Assche for their intelligible code. 
   
   This file contains the classes used to extend a trail core backward 
   inside the kernel. (see Appendix C.1). They are used to define an 
   iterator as defined in extensionsIterator.h
*/

#include <map>
#include "sbox.h"
#include "state.h"
#include "types.h"
#include "trailCore.h"
#include "extensionsIterator.h"
#include "backwardExtension.h"

typedef GenericExtensionIterator< 
        class BackwardInKernelExtensionPreparation, 
        class InKernelTryteColumns,
        class BackwardInKernelExtensionCache,
        class CostFunctionBackwardInKernelExtension,
        class BackwardInKernelExtension> 
        BackwardInKernelExtensionIterator;

/** For a backward extension (A, B) with B in the kernel used to extend
  * a trail core (C, ...), the state B is chosen box-column by box-column.
  * This class is used to store all the possible values that a particular
  * box-column can take. 
  */ 
class InKernelTryteColumns
{       
    private: 
        /** The x-coordinate of the box-column, 0 ≤ xTryte < 3. */
        const int xTryte; 
        /** The z-coordinate of the box-column, 0 ≤ z < 27. */
        const int z; 
        /** Vector that stores all the possible in-kernel values that the 
          * box-column can take (to be compatible throw the map Subtrytes 
          * with another implicit box-column). The elements of that vector 
          * are sorted in ascending order of cost.(see the class
          * TryteColumnSTCompatible in sbox.h)
          */
        vector<TryteColumnSTCompatible> possibleValues;
        /** True if there is no valid value for the box-column. */ 
        bool empty;
        /** The current value of the box-column is given by possibleValues[index] */ 
        int index;
        /** Sbox used to know the box-columns that are in the kernel and 
          * compatible through Subtrytes with another box-column.
          */ 
        Sbox sbox;
        /** Auxiliary variable for the constructor and the method setTryteColumn */ 
        vector<pair<int, int>> tryteAfterST;
    public:
        /** The constructor. It does not initilize possibleValues */
        InKernelTryteColumns(int xTryte, int z):xTryte(xTryte), 
        z(z), index(0){}
        /** @return true if it do not exist a box-column in the kernel and 
          *         compatible with the particular tryte column after 
          *         Subtrytes, false otherwise.
          */
        bool isEmpty();
        /** It sets  the first box-column value.
          * @param cache Not useful here.
          */ 
        void setFirstValue(const BackwardInKernelExtensionCache& cache);
        /** It sets (if possible) the next value of the box-column. 
          * @return true another value has been set, false otherwise.
          */
        bool setNextValue(const BackwardInKernelExtensionCache& cache);  
        /** @return The weight of the current box-column. */
        Weight getWeight() const; 
        /** @return The Hamming weight of the current box-column. */
        int getHammingWeight() const;
        /* It initilizes the vector possibleValues.
         * @param state The state after the map Subtrytes. The box-column must
         *              be compatible with the box-column of @state. 
         */ 
        void setValues(const TroikaState& state);
        /** It sets the current box-column to @state. 
          * @param state The TroikaState to set the box-column to. 
          */
        void setTryteColumn(TroikaState &state) const; 
        friend ostream & operator << (ostream& fout, const InKernelTryteColumns& col);
};

/** Given the activity-pattern of the state C, if first verifies if all the 
  * active box-columns of C have 2 or 3 active trytes. If not, there do not exist 
  * an extension in the kernel. If it may exist a backward in kernel extension, 
  * the method prepareExtension is called to initialize the attributes needed 
  * for the extension. 
  */ 
class BackwardInKernelExtensionPreparation
{
    public:
        long double maxWeightExtension;
        bool possible;
        /** Attribute for the initialization of the vector partsList of the 
          * BackwardInKernelExtensionIterator.
          */ 
        vector<InKernelTryteColumns> tryteColumnsAtB;
        /** Attribute that is going to be used by CostFunctionBackwardInKernelExtension. */
        vector<Weight> toAddToTheCost;
        /** Attribute that is going to be used by CostFunctionBackwardInKernelExtension. */
        int toSubtractFromTheCost;
    public: 
        BackwardInKernelExtensionPreparation(long double aMaxWeightExtension, 
                                             const ActiveState& activeC);
        void prepareExtension(const TroikaState &stateC);
        bool couldBeExtended() const; 
        vector<InKernelTryteColumns>& getPartsList();
        friend ostream & operator << (ostream& fout, const BackwardInKernelExtensionPreparation& prep);
};

/** Auxiliary class for the class BackwardInKernelExtensionIterator.
  * It stores information to compute the cost of the extension that is 
  * being chosen.
  */ 
class BackwardInKernelExtensionCache
{
    public: 
        Weight wBC; 
        unsigned int hammingWeightA; 
    public: 
        BackwardInKernelExtensionCache(BackwardInKernelExtensionPreparation& prep, const TroikaState& stateC);
        void push(const InKernelTryteColumns& tryteColumn);
        void pop(const InKernelTryteColumns& tryteColumn);
};

/** Class used by the BackwardInKenrelExtensionIterator to test if the 
  * extension that is being formed or is already formed has a too high weight. 
  */ 
class CostFunctionBackwardInKernelExtension
{
    public: 
        long double maxWeightExtension; 
        /** Attributes used to compute the cost of a node, as explained in 
          * Appendix C.1
          */
        vector<Weight>& toAddToTheCost; 
        int toSubtractFromTheCost;
    public:
        CostFunctionBackwardInKernelExtension(BackwardInKernelExtensionPreparation &prep, const TroikaState& stateC); 
        bool tooHighCost(const BackwardInKernelExtensionCache& cache, int indCurPart) const;
};

/** The output representation for a backward extension. */
class BackwardInKernelExtension: public BackwardExtension 
{
    public:
        BackwardInKernelExtension(BackwardInKernelExtensionPreparation& prep, const TroikaState& stateC) {} 
        /** This method is called when all the parts of the extension are chosen.
          * It updates the attributes of the BackwardExtension using information 
          * stored in @trytesInfoAtB and @cache. 
          */ 
        void set(const vector<InKernelTryteColumns> tryteColumnsAtB,
                 const BackwardInKernelExtensionCache& cache, 
                 const CostFunctionBackwardInKernelExtension& costF);


};

template<>                     
bool BackwardInKernelExtensionIterator::next(); 
#endif
