/** mixedStateIterator.h */ 
#ifndef MIXED_STATE_ITERATOR_H
#define MIXED_STATE_ITERATOR_H 

/*
This code uses functions and ideas from KeccakTools
(https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

The classes of this file define the tree used to generate the mixed state of the: 
 1) 2-round in kernel trail cores (see Section 5.1)
 2) 2-round out kernel trail cores (see Section 5.2)
*/

#include <stack>
#include "state.h"
#include "traversal.h"
#include "troikaStateIterator.h"
#include "trailCore.h"
#include "forwardExtension.h"

class ActiveTrailCoreCache;
class MixedTrailCoreCache;
class TroikaColumns;
class BareState;
class TroikaStateIterator;
class TwoRoundTrailCoreCostFunction;

typedef GenericTreeIterator<class ActiveTrits,
                            class ActiveTritsSet,
                            class ActiveTrailCoreCache,
                            class TwoRoundTrailCore, 
                            class TwoRoundTrailCoreCostFunction> 
        K_TrailCore_Iterator;

typedef GenericTreeIterator<class Column,
                            class ColumnsSet,
                            class BareStateCache,
                            class BareState, 
                            class TwoRoundTrailCoreCostBoundFunction> 
        BareStateIterator;

typedef GenericTreeIterator<class ActiveTrits,
                            class ActiveTritsSet,
                            class MixedTrailCoreCache,
                            class TwoRoundTrailCore, 
                            class TwoRoundTrailCoreCostFunction> 
        N_TrailCore_Iterator;

/** This class represents the unit of the K_TrailCore_Iterator and N_TrailCore_Iterator. 
  * It contains the column position and y-coordinate(s) of the active 
  * coordinate(s) to add to the mixed state at the input of theta. 
  */
class ActiveTrits: public ColumnPosition
{
    public: 
        /** The 3 least significant bits of isYiActive indicate the 
          * y-coordinate(s) of the active coordinate(s) to add to the mixed
          * state. For 0 ≤ i < 3, if (isYiActive & (1 << i)) != 0, an active
          * trit of coordinates (x, i, z) is added to the mixed state.
          */
        UINT8 isYiActive; 
    public: 
        ActiveTrits() : ColumnPosition(), isYiActive(0) {} 
        /** It increments the position of the column, respecting the lexicographic 
          * order [z, x].
          * @param firstUnit     The first unit of the unit-list 
          * @param optimization  True if the iterator is of type K_TrailCore_Iterator. 
          *                      In that case, coordinates of the unit are incremented 
          *                      so that the unit-list can be canonical.
          * @return false if the coordinates of the column cannot be incremented, 
          *         true otherwise.
          */ 
        bool incrementCoordinates(const ActiveTrits& firstUnit, bool optimization);
        friend ostream & operator << (ostream& fout, const ActiveTrits& trits);
}; 

/** This class represents the set of ActiveTrits and defines the order 
  * relation among them. This set is used by a K_TrailCore_Iterator to 
  * to generate a mixed stated which indicate the active trits positions of 
  * in-kernel Troika states (see Section 5.1). It is also used by a
  * N_TrailCore_Iterator add active trits to the unaffected columns of a mixed 
  * state (see Section 5.2). 
  */
class ActiveTritsSet
{
    public: 
        /** It indicates for each column of the mixed state the first value of
          * isYiActive that can be used by getFirstChildUnit to add active trits 
          * to the column. If this value is 0, that means that active trits 
          * cannot be added to the column. For 0 ≤ x < 9, 0 ≤ z < 27,
          * firstActiveTritsAllowed[x + 9*z] concerns the column of coordinates
          * x and z.
          */
        vector<UINT8> firstActiveTritsAllowed; 
    public: 
        /** The default constructor, used for in-kernel states. */
        ActiveTritsSet():firstActiveTritsAllowed(COLUMNS*SLICES, 0x3){} 
        /** The constructor, for out-kernel states, with firstActiveTritsAllowed
          * specified.
          * @param aFirstActiveTritsAllowed First values of isYiActive that can
          *                                 be set to complete a column. 
          */
        ActiveTritsSet(const vector <UINT8>& aFirstActiveTritsAllowed): 
            firstActiveTritsAllowed(aFirstActiveTritsAllowed){}
        /** It returns an ActiveTrits in the first available position. 
          * @param unitList The list of units used to know the last unit. 
          * @param cache    Not used here.
          * @return The first available ActiveTrits.
          */
        ActiveTrits getFirstChildUnit(
                            const vector<ActiveTrits>& unitList, 
                            const ActiveTrailCoreCache& cache) const; 
        /** It iterates the current unit.
          * @param unitList Not used here. 
          * @param current  The current unit.  
          * @param cache    Not used here.
          */ 
        void iterateUnit(const vector <ActiveTrits>& unitList,
                         ActiveTrits& current,
                         const ActiveTrailCoreCache& cache) const;
        /** It checks if a list of ActiveTrits is z-canonical
          * (only for K_TrailCore_Iterator).
          */
        bool isCanonical(const vector<ActiveTrits>& unitList,
                         ActiveTrailCoreCache& cache) const;
    private: 
        /** It compares two given ActiveTrits, using the lexicographic
          * order [z, x, isYiActive].
          * @param first   The first ActiveTrits.
          * @param second  The second ActiveTrits.
          * @return 0 if the two ActiveTrits are equal, 
          *         1 if the first is smaller,
          *         2 if the second is smaller.
          */ 
        unsigned int compare(const ActiveTrits& first, 
                             const ActiveTrits& second) const;
}; 

/** This class is used during the traversal of a K_TrailCore_Iterator or a 
  * N_TrailCore_Iterator to represent the active trits' positions of the 
  * 2-round trail core A --Lambda--> B that is beeing formed. 
  * It contains all the information needed about the active 
  * state A and the active state B to compute the cost of the trail core.
  */
class ActiveTrailCoreCache
{
    public: 
        /** The positions of the active trits of A, completed during the 
          * traversal of the K_TrailCore_Iterator or N_TrailCore_Iterator.  
          */
        ActiveState activeA;
         /** The positions of the active trits of A, completed during the 
          * traversal of the N_TrailCore_Iterator.  
          */
        ActiveState activeB; 
        /** The stack for the minimum reverse weight of A. */ 
        stack <unsigned int> wA;  
        /** The stack for the minimum direct weight of B. */
        stack <unsigned int> wB; 
        /** It indicates whether the last push was a dummy push or not. */
        bool dummy; 
        /** True if the cache is used for a K_TrailCore_Iterator, false if the
          * cache is used for an N_TrailCore_Iterator. 
          */ 
        bool kernel;
        /** The position of the unaffected column of parity zero of the mixed 
          * state (see troikaStateIterator.h ). Their attributes columnValues
          * must be initialized later.
          */ 
        vector<TroikaColumns> inKernelColumns; 
    public:
        /** A constructor. It is used to initialize the cache of a 
          * K_TrailCore_Iterator.
          */
        ActiveTrailCoreCache(); 
        /** It performs a dummy push to the cache. */
        void pushDummy(); 
        /** It pushes to the cache the trits used to complete a column.
          * It updates the attributes dummy, activeA, activeB,  
          * wA and wB.
          * @param activeTrits The trits to be pushed.
          */
        void push(const ActiveTrits& activeTrits);
        /** It pops the highest unit from the cache. It updates activeA, 
          * activeB, wA and wB.
	      * @param activeTrits The activeTrits to pop.
	      */
        void pop(const ActiveTrits& activeTrits);
        friend ostream & operator << (ostream& fout, const ActiveTrailCoreCache& cache);
};

/** This class is used during the traversal of a N_TrailCore_Iterator to
  * represent the active trits' positions of the 
  * 2-round trail core A --Lambda--> B that is beeing formed. 
  * It contains all the information needed about the active 
  * state A and the active state B to compute the cost of the trail core.
  */
class MixedTrailCoreCache : public ActiveTrailCoreCache
{
    public:
        /** The TroikaColumns chosen during the traversal of a 
          * BareStateIterator. The affected columns are fixed, their 
          * attributes columnValues are already initialized. The 
          * unaffected columns can by completed during the traversal of 
          * the N_TrailCore_Iterator, their attributes columnValues
          * must be initialized later. 
          */
        vector<TroikaColumns> outKernelComponentColumns;
        /** The number of components of the state at the input of theta. */ 
        unsigned int nrComponents;
    public:
        MixedTrailCoreCache(const BareState& parityBareState);
        friend ostream & operator << (ostream& fout, const MixedTrailCoreCache& cache);
};

/** The output representation of a TwoRoundTrailCore (A, B) generated during the
  * traversal of a N_TrailCore_Iterator or a K_TrailCore_Iterator. 
  */ 
class TwoRoundTrailCore
{
    public:
        /** The minimum reverse weight of the state A. */ 
        unsigned int wA; 
        /** The minimum direct weight of the state B. */ 
        unsigned int wB; 
        TroikaStateIterator statesB;
        /** The trit activity pattern of the state A. */ 
        ActiveState activeA;
        /** The trit activity pattern of the state B. */ 
        ActiveState activeB;
    public: 
        TwoRoundTrailCore(): wA(0), wB(0){};
        /** This methods outputs the 2-round trail core to save it in, e.g., a file.
	      * @param fout The stream to save the trail core to.
	      */
        void save(ostream& fout);
        /** This methods sets the trail. It is used by a K_TrailCore_Iterator.
       	  * @param unitList The list of activeTrits representing the trail.
          * @param cache The cache representation of the trail.
          * @param costF Not used here. 
          * @param aMaxCost Not used here.
	      */
        void set(const vector<ActiveTrits>& unitList,
                 const ActiveTrailCoreCache& cache,
                 const TwoRoundTrailCoreCostFunction& costF, 
                 unsigned int aMaxCost); 
        /** This methods sets the trail.
	      * @param unitList The list of activeTrits representing the trail.
	      * @param cache The cache representation of the trail.
          * @param costF  Not used here.
          * @param aMaxCost Not used here.
 	      */
        void set(const vector<ActiveTrits>& unitList,
                 const MixedTrailCoreCache& cache, 
                 const TwoRoundTrailCoreCostFunction& costF, 
                 unsigned int aMaxCost); 
}; 

/** This class is used to compute a bound on the cost of a 2-round trail core
  * generated during the traversal of a K_TrailCore_Iterator or a 
  * N_TrailCore_Iterator and of its children.
  * The cost of a 2-round trail core (A, B) is alpha * wMinRev(A) + beta * wMinDir(b)
  */
class TwoRoundTrailCoreCostFunction {

    public:
        unsigned int alpha;
        unsigned int beta;
    public:
        TwoRoundTrailCoreCostFunction()
            : alpha(1), beta(1){}
        TwoRoundTrailCoreCostFunction(unsigned int aAlpha, unsigned int aBeta)
            : alpha(aAlpha), beta(aBeta) {}
	    /** It returns the bound on the cost of the trail core and 
          * its children.
          */ 
        unsigned int getCost(const vector<ActiveTrits>& unitList,
                             const ActiveTrailCoreCache& cache) const;
};

/** This function save in @fout all the 2-round trail cores (A, B)
  * of cost @alpha * wMinRev(A) + @beta * wMinDir(b) below @aMaxCost.
  */ 
void traverse_K_TrailCoresTree(unsigned int aMaxCost,
                               unsigned int alpha,
                               unsigned int beta, 
                               ostream &fout); 
#endif
