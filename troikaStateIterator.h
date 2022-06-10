/* troikaStateIterator.h */
#ifndef TROIKA_STATE_ITERATOR_H
#define TROIKA_STATE_ITERATOR_H

/*
This code uses functions and ideas from KeccakTools
(https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

The classes of this file are used to form, from a mixed state of Troika 
(see Section 5), all the Troika states that respect some conditions 
(eg conditions on the parity of the state, trits that are necessarily active ...).
These states are chosen column by column. 
*/

#include "state.h"

class  MixedTrailCoreCache;
class ActiveTrailCoreCache; 

#define NR_IN_KERNEL_COLUMN 9

/** Class used to choose, from a list of possible values, the value of a column. */
class TroikaColumns : public ColumnPosition 
{
    public:
        /** It indicates if the column is affected or not. */ 
        bool isAffected; 
        /** It indicates the parity of the column. */ 
        unsigned int parity;
        /** Index of the component of the parity bare state (see Section 5.5) 
          * of the column. If indexOutKernelComponent = -1, that means that the
          * column does not belong to a component : it is an unaffected column 
          * of parity 0.
          */  
        int indexOutKernelComponent; 
        /** An array<int, 3> is seen as a column value : it stores the 3 
          * values of the 3 trits of a column. The vector
          * possibleColumnValues stores all the values that can be taken 
          * by the column. 
          */
        vector<array<int, 3>> columnValues;
        /** possibleColumnValues[indexValue] is the current column value. 
          */
        unsigned int indexValue; 
    public:
        /** The constructor. It does not initialize the attribute columnValues. 
          * This attribute must be initialized with the functions setValues 
          * and addValues.
          */ 
        TroikaColumns(unsigned int aX, unsigned int z, bool affected = false, 
                      unsigned int parity = 0, int componentIndex = -1); 
        /** It sets (if possible) the next column value.
          * @return false if all the possible values have already been chosen,
          *         true otherwise.
          */
        bool next();
        /** It returns the current trit value of the column of coordinate y.
          * @param y The coordinate, 0 ≤ y < 3, of the trit value to get.  
          * @return  The trit value of coordinate @y of the current column value.  
          */
        unsigned int getTrit(unsigned int y) const;
        /** It uses the vector @values to initialize the attribute columnValues. */ 
        void setValues(const vector<array<int, 3>>& values);
        /** It adds to the vector columnValues the array @value*/ 
        void addValues(const array<int, 3>& value); 
        friend ostream & operator << (ostream &fout, const TroikaColumns &aColumn);        
}; 

/** This class is used to : 
  * 1) form all the in-kernel states that satisfy conditions relating to 
  *    coordinates necessarily active or necessarily passive (see 
  *    Appendix C.3 Forward in kernel extension and forwardInKernelExtension.h) 
  * 2) form all the in-kernel states compatible with a given mixed state 
  *    (see Section 5.1 Generating |K|-trail cores).
  * 3) form all the Troika states compatible with a given mixed states 
  *    and a given parity plane.
  */
class TroikaStateIterator
{
    public: 
        /** A list of the TroikaColumns associated to a mixed state. 
          * This list is used to form all the Troika states 
          * by taking all combinations of column values. 
          */
        vector<TroikaColumns> columnsList; 
        /** This attribute is used during the tree traversal. It is the index, 
          * in the vector columnsList, of the column whose value must be chosen.
          */
        int indexColumn;
        /** Variable used to store the current state. */ 
        TroikaState state; 
        /** It indicates if the iterator has reached the end. */
        bool end; 
        /** This attribute is only used when the TroikaStateIterator is
          * associated to an out-kernel mixed state. It indicates the number
          * of components of the state (see Section 5.5). 
          */
        unsigned int nrOutKernelComponents;
        /** This attribute is only used when the TroikaStateIterator is
          * associated to an  out-kernel mixed state. During the traversal
          * of the tree, it is incremented (from 0 to 2**nrComponents - 1). 
          * The bits of the integer multiplyComponentBy2 indicates which 
          * component has to be multiplied by 2.
          */
        unsigned int multiplyOutKernelComponent;
        /** Array used to store the 9 in kernel column values. */
        const static array<int, 3> inKernelColumn[9];
        /** The attribute possibleColumnValues is used to choose the
          * unaffected columns of a mixed state.
          * For 0 ≤ sumsTo < 3, for 0 ≤ isYiActive < 8, 
          * possibleColumnValues[sumsTo][isYiActive] contains a vector 
          * of array<int, 3>. The columns / array<int, 3> stored by the vector
          * are such that : 
          * - the parity of the column equals sumsTo
          * - the trit of the column of coordinate y is active iff 
          *   (isYiActive & (1 << y)) != 0.  
          */
        static vector<array<int, 3>> **possibleColumnValues;
    public:
        TroikaStateIterator(): end(true){}
        /** A constructor. It is used to form all the in-kernel states
          * compatible with the mixed state @activeInKernel.
          * @param activeInKernelState The columns of this active state must 
          *                            have 0, 2 or 3 active trits.
          */
        TroikaStateIterator(const ActiveState& activeInKernel);
        /** A constructor called during a forward in kernel extension.
          * It initializes the columnList with valid Troika columns. A column
          * of a Troika state is valid if : 
          * - it is in the kernel, 
          * - it respects the trits necessarily active, indicated by
          *   @mandatoryActiveTrits
          * - it respects the trits necessarily passive, indicated by 
          *   @mandatoryActiveTrits.
          * @param possibleActiveTrits  Variable which indicates the trits  
          *                             necessarily passive of a Troika state. 
          *                             It is the trits of possibleActiveTrits
          *                             equal to 0.
          * @param mandatoryActiveTrits Variable which indicates the trits
          *                             necessarily active of a TroikaState. 
          *                             It is the trits of mandatoryActiveTrits
          *                             equal to 1.
          */
        TroikaStateIterator(const ActiveState& possibleActiveTrits, 
                            const ActiveState& mandatoryActiveTrits);
        /** A constructor. It is used to generate all the in-kernel states 
          * compatible with the current node of a K_TrailCore_Iterator. 
          */
        TroikaStateIterator(const ActiveTrailCoreCache& cache);
        /** A constructor. It is used to generate all the out-kernel states 
          * compatible with the current nod of a N_TrailCore_Iterator.
          */
        TroikaStateIterator(const MixedTrailCoreCache& cache);
        /** @return True if there are no more states to consider. */
        bool isEnd() const;
        void first(); 
        /** It moves the iterator to the next state. */
        void operator++();
        /** It gives a constant reference to the current state. */
        const TroikaState& operator*();
        friend ostream & operator << (ostream &fout, const TroikaStateIterator &it);     
    private:  
        bool toChild();
        bool toSibling(); 
        bool toParent();
}; 

#endif
