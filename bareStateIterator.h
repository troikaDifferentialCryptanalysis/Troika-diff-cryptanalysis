/** bareStateIterator.h */
#ifndef BARE_STATE_ITERATOR_H
#define BARE_STATE_ITERATOR_H

/*
This code uses functions and ideas from KeccakTools
(https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

The classes of this file define the tree used to generate the 2-round trail 
cores (A --> rho --> theta --> B) outside the kernel such that rho(A) is a
parity-bare (see section 5.3).
*/


#include <vector>
#include <algorithm>
#include <iostream>
#include <list>
#include <queue>
#include <map>
#include <assert.h>
#include "types.h"
#include "state.h"
#include "traversal.h"
#include "troikaStateIterator.h"

class BareStateCache;

/** If there is already in the unit-list another column of same 
  * coordinates, the column is considered entangled. Otherwise, the column's 
  * type is noEntanglement.
  *
  * There are 2 kinds of authorized entanglements : 
  * - an affected column of parity 0 overlapping a column with only one
  *   active trit of y-coordinate = 0, 
  * - a column with only one active trit of y-coordinate = 0 overlapping
  *   an affected column of parity 0.
  * In these cases, the column's type can be noNeedToIterate, 
  * toIterate or iterated. If the entangled column connects two supra-units that
  * already belonged to the same component of the supra-unit graph 
  * (see section 5.5), the column's type is noNeedToIterate. Otherwise, the
  * column's type is first toIterate, then the method iterateEntanglement is
  * called and the column's type becomes iterated.
  */
enum Entanglement {noEntanglement, 
                   noNeedToIterate,
                   toIterate, 
                   iterated};

/** The unit of the tree. See paragraph "From supra-units to units" of 
  * Section 5.3. 
  */
class Column
{
    public:
        /* The d-coordinate of the supra-unit of the column. */ 
        int dSupraUnit; 
        /* The z-coordinate of the supra-unit of the column. */ 
        int zSupraUnit; // It has to be an int (not an unsigned int) for the method zTranslateSupraUnit  
        /** The index of the unit in the supra-unit (see the end of the 
          * paragraph "From supra-units to units"). The unit that starts the 
          * supra-unit has rank = 0. If the rank is even, then then the unit 
          * is affected. 
          */ 
        unsigned int rank;
        /** True if the unit ends a supra-unit, false otherwise. */ 
        bool ending;
        /** Index that indicates the value of the column before theta.
          * The values are defined up to the multiplication by the scalar 
          * @scalar.
          * An affected unit can take 9 values. These values are stored 
          * in the array ZERO_PARITY_COLUMN. If the column is affected, its 
          * value is mod 3 : 
          * - before theta : ZERO_PARITY_COLUMN[indexValue] * scalar 
          * - after theta  : (ZERO_PARITY_COLUMN[indexValue] + 1 ) * scalar.
          * Thus, two affected columns with the same index value and different 
          * scalar have their active trits at the same position.
          * If the column is unaffected, it has before and after theta 
          * only one active trit of y-coordinate given by @indexValue and of 
          * value given by @scalar. 
          */ 
        unsigned int indexValue;
        /** The x-coordinate of the column. */ 
        unsigned int x; 
        /** The z-coordinate of the column. */
        unsigned int z; 
        /** see Entanglement. */ 
        enum Entanglement entanglementType;
        /** Used to know the value of the column (see index value). */ 
        unsigned int scalar; 
    public:
        Column(); 
        void incrementRank(); 
        bool incrementIndexValue();
        /** It increments the supra-unit position, respecting the lexicographic 
          * order [d, z, indexValue]. For the canonicity, the coordinates of 
          * the first supra-unit must satisfy dSupraUnit = 0 and
          * 0 ≤ zSupraUnit < 9.
          */ 
        bool incrementSupraUnitPosition(const vector<Column>& unitList); 
        void changeScalar(); 
        unsigned int affectedBy() const; 
        unsigned int parity() const;
        /** It returns the z-coordinate of the smallest column for the 
          * lexicographic order [d, z] that has the same x-coordinate than 
          * the column of (d, z)-coordinates (dSupraUnit, zSupraUnit).
          * The coordinates returned must be between 0 and 8. 
          */ 
        unsigned int zSupraUnitCanonical() const;
        /** Auxiliary method for the method isCanonical of the ColumnsSet.
          * It updates the attributes dSupraUnit and zSupraUnit with the 
          * (d, z)-coordinates of the column that has its z-coordinate equal 
          * to (zSupraUnit + dz) and its x-coordinate equal to the x-coordinate 
          * of the column of (d, z)-coordinates (dSupraUnit, zSupraUnit).
          */ 
        void zTranslateSupraUnit(int dz);
        /** An ordering operator. It is the lexicographic order 
          * [dSupraUnit, zSupraUnit, rank, (ending = false) < (ending = true), indexValue]. 
          * It is not a total ordering. Indeed, columns are equivalent up to 
          * a multiplication by the scalar 2. 
          */ 
        bool operator < (const Column& other) const;
        friend ostream & operator << (ostream& fout, const Column& col);
    private:
        void setCoordinates();
}; 

/**  This class represents the set of column assignments and defines the 
  *  order relation among them.
  */
class ColumnsSet
{
    public: 
   	ColumnsSet(){ };  
	Column getFirstChildUnit(vector<Column>& unitList,
                             BareStateCache& cache) const;
	/** This method iterates the current unit with respect to the order 
      * relation [dSupraUnit, zSupraUnit, rank, ending, indexValue]
      * and restrictions on its type (scalar, entanglement). It iterates if 
      * needed an entanglement by multiplying by 2 a component. 
	  */
	void iterateUnit(vector<Column>& unitList, Column& current, BareStateCache& cache) const;	
	/** This method checks if a list of column assignments is z-canonical.
      * The order relation is not a total order relation.  
      * Two bare states that have the same supra-unit up to a scalar are 
      * considered equal. 
      */
	bool isCanonical(const vector<Column>& unitList, const BareStateCache& cache) const;
private:
	/** This method checks if a given column overlaps any other column of
      * the unit-list. If it is true and the overlaps is not allowed, the column
      * cannot be added to the unit list. If it is false, the column
      * can have either the attribute toIterate, noNeedToIterate or iterated.
      * The attribute of the column is "noNeedToIterate"" if adding the column 
      * do not "merge" two components of the bare state (see section 5.5 for 
      * the definition of a component). Otherwise, the column is first
      * "toIterate". It becomes "iterated" when the method iterateEntanglement
      * is called.
      *
      * @param current  A column. We want to know if that column can be added
      *                 to the unit list. 
      * @param cache    A reference to the cache representing the trail core. 
      *                 It is used to know if @current overlaps a unit already
      *                 in the unit list.
      * @return true if the column can be added to the unit-list, false
      *         otherwise.
	  */
	bool checkColumnOverlapping(const vector<Column>& unitList, Column& current,
                                const BareStateCache& cache) const;
};
/** Class used to store a pair of coordinates. */
class TritPositionAtAAndB
{
    public: 
        /** coordinate at A. */ 
        TritPosition pA; 
        /** coordinate at B. */ 
        TritPosition pB; 
    public: 
        TritPositionAtAAndB(const Column& col, unsigned int y)
        : pB(col.x, y, col.z) { pA = pB.getInvSRSL(); }
        TritPositionAtAAndB(unsigned int x, unsigned int y, unsigned int z)
        : pB(x, y, z) { pA = pB.getInvSRSL(); }

};

/** We denote be (A, B) the 2-round trail cores generated by the tree traversal.
  * This class stores information about a particular trit of a unit at A and B.
  */
class TritAtAAndB : public TritPositionAtAAndB
{
    public:
        /* The value of the trit at A. */ 
        unsigned int vA; 
        /* The value of the trit at B. */ 
        unsigned int vB;
        /** The d-coordinate of the supra-unit that contains that trit. */ 
        int dSupraUnit; // int (and not unsigned int) 
        /* true if we know that the trit is stable. A trit can turned stable 
         * if new units are added to the unit-list. 
         */ 
        bool stable; 
    public: 
        TritAtAAndB(const Column& col, unsigned int y); 
        friend ostream & operator << (ostream& fout, const TritAtAAndB& trit); 
};


/** This class represents a 2-round trail core A --Lambda--> B, where 
  * rho(A) is a parity-bare state.
  * It contains all the information needed about A, B, the runs of rho(A) and 
  * the cost of the trail core during the traversal of the tree.
  */
class BareStateCache
{
    public: 
        TroikaState  stateA; 
        TroikaState  stateB; 
        unsigned int nrActiveTrytesA;
        unsigned int nrActiveTrytesB;
        /** Trit-activity pattern that indicates some trits of stateA that have already been identified as safe. */ 
        ActiveState  stableTritsA; 
        /** Trit-activity pattern that indicates some trits of stateB that have already been identified as safe. */ 
        ActiveState  stableTritsB; 
        /** Number of active trytes of stableTritsA. */ 
        unsigned int nrStableTrytesA; 
        /** Number of activity trytes of stableTritsB. */ 
        unsigned int nrStableTrytesB;
        /** The parity plane of rho(A). */
        TroikaPlane  parityPlane; 
        /** The theta-effect plane of rho(A). */ 
        TroikaPlane  thetaEffect; 
        vector<TritAtAAndB> possibleUnstableTrits; 
        /** It indicates whether the last push was a dummy push or not. */
        bool dummy; 
        /** startingSupraUnit[i] contains the index of the unit that starts the 
          * i-th supra unit of the unit-list. (The smallest index is 0).
          */
        vector<unsigned int> startSupraUnit; 
        /** This array is used during the traversal of the tree to know the 
          * index of the supra-unit of the unit. 
          * An entangled column can belong to 2 supra-units. 
          * In that case, the array stores the index of the supra-unit that
          * was added first. For 0 ≤ x < 9,  0 ≤ z < 27,
          * supraUnitIndexes[x][z] = -1 means that none of the column has 
          * coordinates x and z. 
          */
        int supraUnitIndexes[COLUMNS][SLICES];
        /** Index of the last supra-unit of the unitList. */ 
        int indexLastSupraUnit;
        /** Index of the last unit added to the unit-list. */ 
        int indexLastUnit;
        /** Graph used to organize the components of rho(A). See Section 5.5. */ 
        vector<list<int>> neighborsSupraUnit;
    public:
        BareStateCache();  
        void pushDummy(); 
        void push(const Column& col); 
        void pop(const Column& col);
         /** An entanglement can lead to : 
          * - a column of parity 1 ≤ p < 3 and affected by p, 
          * - a column of parity 1 ≤ p < 3 and affected by (2p)%3;
          * This entanglement is formed by a unit_1 which belongs to a
          * component_1 of rho(A). When the entangled column unit_2 
          * overlaps unit_1, the component_1 of the run graph might be merged 
          * with another component of rho(A). In that case, the attribute 
          * Entanglement of unit_2 is "toIterate". This method multiplies 
          * by 2 the columns of the component component_1. After that, the 
          * attribute of unit_2 becomes "iterated".
          * @param unitList  The unit list containing the columns. 
          * @param aColumn   The column whose attributes Entanglement is
          *                  "toIterate".
          */
        void iterateEntanglement(vector<Column> &unitList, Column& current);  
        /** It indicates if the last supra-unit and the supra-unit of index 
          * @index are in the same component.
          */
        bool lastSupraUnitInTheComponentOf(unsigned int index) const;
        /** It returns a vector of the form <C1, C2, ..., Cn> where Ci 
          * is a vector that stores all the supra-unit indexes of a 
          * component of rho(A).
          */ 
        vector<vector<unsigned int>> getUnitIndexesOfAllTheComponents() const; 
    private:
        /** It gets the indexes of a component and indicates in the vector
         * @supraUnitFound that these indexes have been found. 
          * @param supraUnitIndexOfTheComponent An index of a component of 
          *                                     rho(A), between 0 and
          *                                     neighborsSupraUnit.size() - 1
          * @param supraUnitFound  A vector whose size is equal to the number 
          *                        of components of rho(A). 
          * @return A vector whose elements are the indexes of the supra-unit 
          *         that are in the same component as the supra-unit of index 
          *         @supraUnitIndexOfTheComponent. 
          */ 
       vector<unsigned int> getUnitIndexesOfAComponent(
                                unsigned int supraUnitIndexOfTheComponent,
                                vector<bool>& supraUnitFound) const;
       /** Auxiliary function for the method push and pop that updates
         * the attributes nrActiveTrytesA, stateA, nrActiveTrytesB, stateB, 
         * nrStableTrytesA, nrStableTrytesB, stableTritsA, stableTritsB and 
         * possibleUnstableTrits. 
         */ 
        void pushOrPopTritAtAAndB(bool push, const TritAtAAndB& trit);
        void addSupraUnit();
        void removeLastSupraUnit();
        /** It adds to the graph used to organize the components of rho(A) an 
          * edge between the index of the current supra-unit and @index.
          */ 
        void addSupraUnitNeighbor(unsigned int index);
        /** It removes from the graph used to organize the components of rho(A)
          * the edge between the index of the current supra-unit and @index.
          */ 
        void removeSupraUnitNeighbor(unsigned int index);
};

/** Class used to compute a lower bound on the costs 
  * alpha * wMinRev(A) + beta * wMinDir(B) of a nodes (A, B) and its
  * its descendants. 
  */ 
class TwoRoundTrailCoreCostBoundFunction
{
    public: 
        unsigned int alpha; 
        unsigned int beta;
    public: 
        TwoRoundTrailCoreCostBoundFunction(unsigned int aAlpha, unsigned int aBeta)
            : alpha(aAlpha), beta(aBeta){};
        /** It returns a lower bound on the costs of a unit-list and its 
          * descendants computed using Algorithm 1 of Appendix A.
          */ 
        unsigned int getCost(const vector<Column>& unitList, 
                             const BareStateCache& cache) const;
    private:
        bool isTritStillUnstable(const TritAtAAndB& trit,
                                 const BareStateCache& cache,
                                 const vector<Column>& unitList) const; 
        unsigned int getContributionUnstableTrit(
                                 const TritAtAAndB& trit, 
                                 ActiveState& possibleActiveTritsA, 
                                 ActiveState& possibleActiveTritsB) const; 
};

/** The output representation of a trail core (A, B). */ 
class BareState
{
    public: 
        TroikaState stateA; 
        TroikaState stateB;
        /** The minimum reverse weight of A. */ 
        unsigned int wA; 
        /** The minimum direct weight of B. */ 
        unsigned int wB;
        /** */ 
        vector<TroikaColumns> outKernelComponentsColumns;  
        /** Number of components of the parity bare state rho(A). */ 
        unsigned int nrComponents;
        /** It indicates for each column of the parity bare state the first
          * which trits can be activated by a N_TrailCore_Iterator
          * (see mixedStateIterator.h). 
          */
        vector<UINT8> firstActiveTritsAllowed;
        /** True if the unit-list is complete and if the cost of the node is
          * not too high. 
          */ 
        bool valid; 
    public : 
        BareState(){}
        /** It indicates if the current node is a valid trail core. If it is 
          * the case, the attributes of the current BareState are updated. 
          */ 
        void set(const vector<Column>& unitList, 
                 const BareStateCache& cache, 
                 const TwoRoundTrailCoreCostBoundFunction& costF, 
                 unsigned int aMaxCost); 
        friend ostream & operator << (ostream& fout, const BareState& bareState); 
};

#endif
