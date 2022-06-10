/* KK_trailCores.h */ 
#ifndef KK_TRAILCORES_H
#define KK_TRAILCORES_H

/* 
   This code uses functions and ideas from KeccakTools
   (https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
   and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
   We thank Silvia Mella and Gilles Van Assche for their intelligible code. 
  
   The classes of this file are used to generate all 3-round trail cores of 
   parity profile |K|K| up to some weight (see Appendix B).
   The 3-round trail cores are denoted by (A, B, C, D). The first step is to 
   generate the trit-activity patterns of the states C and D with a tree traversal.
   
*/

#include <fstream>
#include <iostream> 
#include <iterator>
#include <ostream>
#include <sstream>
#include <map>
#include <set>
#include "state.h"
#include "traversal.h"
#include "sbox.h"
#include "troikaStateIterator.h"
#include "trailCore.h"
#include "backwardInKernelExtension.h"

/** Iterator to generate the trit-activity patterns of the states C and D.
  * For this iterator, the notion of canonicity is not used.
  */
 typedef GenericTreeIterator<class ActiveTritAtCAndD, 
                            class ActiveTritsAtCAndDSet, 
                            class ActiveStatesCAndDCache, 
                            class ActiveStatesCAndD, 
                            class KK_TrailCoreCostFunction> 
        ActiveStatesCAndDIterator; 
 
 /* We say that a unit-list is valid if all the units of this list have a
  * neighbor by constraint on C in the list and a neighbor by constraint on D 
  * in the list (see Appendix B).
  *
  * The unit-list of this iterator can be decomposed into sub-lists of the form  
  * L = L1 || L2 || ... || L{n - 1} || Ln, where for 1 ≤ i < n, 
  * the list L1 || ... || Li is valid. The list Ln can be valid or not. 
  * The units that start a sub-list are called "startingTrit".  
  * The type of the unit indicates if the unit : 
  * - is on the same tryte column at C than its father, 
  * - is on the same column at D than its father,
  * - is on the same tryte column at C than the trit that begins the last 
  *   sub-list Ln (not valid yet) 
  * - is on the same column at D than the trit that begins the last sub-list 
  *   Ln (not valid yet). 
  * - begins a new sub-list of the unit-list.
  * These 5 cases are described in Algorithm 2 of Appendix B.1. 
  *
  */
enum Type {onSameTryteColumnAtC,
           onSameColumnAtD,
           onTheSameTryteColumnAtCAsTheLastStartingTrit,
           onTheSameColumnAtDAsTheLastStartingTrit,
           startingTrit};

/** This class represents a unit of the iterator ActiveStatesCAndDIterator. 
  * It represents an active trit, both from the point of view of the states 
  * C and D.
  */
class ActiveTritAtCAndD
{
    public:
        /** The type of the trit. See Type. */ 
        enum Type type; 
        /** The trit's position at D. */ 
        TritPosition pD;  
        /** The trit's position at C. */ 
        TritPosition pC;
        /** The position of the unit (at C or at D, depending on the type of
          * the unit) set by the method getFirstChild. The method iterateUnit
          * uses this attribute and the attributes yOffset and xOffset to change
          * the position of the unit.
          */ 
        TritPosition firstPosition;
    private:
        /** If the type of the trit is onSameTryteColumnAtC or 
          * onTheSameTryteColumnAtCAsTheLastStartingTrit, then 0 ≤ yOffset < 2
          * and 0 ≤ xOffset < 3.
          * If the type of the trit is onSameColumnAtD or 
          * onTheSameColumnAtDAsTheLastStartingTrit, then xOffset = 0 and 
          * 0 ≤ yOffset < 2.
          * If the trit'type is onTheSameColumnAtDAsTheStartingPositionSubList
          * or onSameColumnAtD, then positionAtD.x = firstPosition + xOffset, 
          *                          positionAtD.y = firstPosition + yOffset
          *                          positionAtD.z = firstPosition.z.
          * If the trit's type is onTheSameTryteColumnAtCAsTheStartingPositionSubList
          * or onSameTryteColumnAtC, then positionAtC.x = firstPosition + xOffset, 
          *                               positionAtC.y = firstPosition + yOffset 
          *                               positionAtC.z = firstPosition.z.
          * If the type of the trit is newStartingPosition, the trit's position
          * at D is given by the attribute firstPosition.
          */ 
        unsigned int yOffset;
        unsigned int xOffset; 
    public:
        ActiveTritAtCAndD();
        bool incrementOffSets();
        /* It uses the attributes type, firstPosition, xOffset and yOffset to 
         * compute pD and pC.
         */ 
        void setCoordinates();
        /** Ordering operator. It is the lexicographic order [x, y, z] on 
          * the trit position at C.
          */ 
        bool operator < (const ActiveTritAtCAndD& other) const;
        friend ostream & operator << (ostream &fout, const ActiveTritAtCAndD& trit);
};

/** This class organizes the tree, as defined in Algorithm 2 of Appendix B.1. */ 
class ActiveTritsAtCAndDSet
{
    public: 
        ActiveTritsAtCAndDSet(){}; 
        ActiveTritAtCAndD getFirstChildUnit(
                            const vector<ActiveTritAtCAndD>& unitList, 
                            const ActiveStatesCAndDCache& cache) const;
        void iterateUnit(const vector<ActiveTritAtCAndD>& unitList, 
                         ActiveTritAtCAndD& current,
                         const ActiveStatesCAndDCache& cache) const;
        /** Here, we do not use the order relation on the units to organize the
          * unit-list.  Therefore the notion  of a "canonical" unit-list does not
          * exist. This method is just used to avoir having the same valid 
          * pattern twice. All the valid patterns are saved and when a new valid 
          * pattern is found, we just check if we already know it or not. 
          */ 
        bool isCanonical(const vector<ActiveTritAtCAndD>& unitList, 
                         const ActiveStatesCAndDCache& cache) const;
    private: 
        enum Type getFirstChildType( 
                          const ActiveTritAtCAndD& parent, 
                          const ActiveStatesCAndDCache& cache) const;
        bool isAValidUnit(const ActiveTritAtCAndD& trit,
                          const vector<ActiveTritAtCAndD>& unitList, 
                          const ActiveStatesCAndDCache& cache) const; 
        ActiveTritAtCAndD getCandidateForFirstChildUnit(
                            const vector<ActiveTritAtCAndD>& unitList, 
                            const ActiveStatesCAndDCache& cache) const;  
};

/** The output representation for the valid activity patterns (C, D). */ 
class ActiveStatesCAndD {
    public: 
        ActiveState activeD;
        ActiveState activeC;
        unsigned int wMinRevC; 
        unsigned int wMinDirD;
    public: 
        ActiveStatesCAndD(): wMinRevC(0), wMinDirD(0){};
        void set(const vector<ActiveTritAtCAndD>& unitList, 
                 const ActiveStatesCAndDCache& cache, 
                 const KK_TrailCoreCostFunction& costF, 
                 unsigned int aMaxCost); 
};

/** This class is used during the traversal of the ActiveStatesCAndDIterator.
  * It contains all the information needed to control the weight during the 
  * tree traversal. (see Algorithm 3 of Appendix B.2). 
  * It also saves all the valid patterns found. To avoir having 
  * the same valid pattern twice, when a new valid pattern is found, we just
  * check if we already know it or not. 
  */
class ActiveStatesCAndDCache
{
    public:
        /** Set that stores the biggest reprensatives of the valid states C
          * and D that were found. It is used to avoid saving states for which 
          * a reprensative is already known.
          */
        set<ActiveState> patternsC; 
        ActiveState stateC; 
        ActiveState stateD;
        /** Variable used to compute a lower bound on the weights of a 
          * 3-round trail core (A, B, C, D), where the activity patterns of
          * C and D are given by the active states stateC and stateD. 
          * In Algorithm 3 of Appendix B.2, this variable is denoted Delta. 
          */ 
        ActiveState possibleTrytesA;
        /** A lower bound on the number of active trytes of the state A. 
          * In Algorithm 3 of Appendix B.2, this variable is denoted min(wBox(Delta)).
          */ 
        int lowestNrActiveTrytesA;
        /** The number of active trytes of stateC. */ 
        unsigned int nrActiveTrytesC; 
        /** The number of active trytes of stateD. */ 
        unsigned int nrActiveTrytesD; 
        /** For  0 ≤ x < 9,  0 ≤ y < 3,  0 ≤ z < 27, 
          * nrTritsOnSameColumnAtD[x][y][z] indicates the number of 
          * neighbor by "constraint on D" that the unit of (x, y, z) coordinates
          * at D has. The unit-list is valid if all of its ActiveTritAtCAndD 
          * have a neighbor by "constraint on D" and a neighbor by 
          * "constraint on C".
          */
        int nrTritsOnSameColumnAtD[COLUMNS][ROWS][SLICES];
        /** For  0 ≤ x < 9,  0 ≤ y < 3,  0 ≤ z < 27, 
          * nrTritsOnSameTryteColumnAtC[x][y][z] indicates the number of 
          * neighbor by "constraint on C" that the trit of (x, y, z) coordinates
          * at C has. 
          */
        int nrTritsOnSameTryteColumnAtC[COLUMNS][ROWS][SLICES]; 
        /** Vector that stores every ActiveTritAtCAndD that begin a new 
          * sub-list (see enum Type to have the description on a sub-list).
          * The first unit of the next sub-list must be bigger (for the 
          * lexicographic order [x, y, z]) than startingPositionAtD.back().
          */
        vector<ActiveTritAtCAndD> startingTrits;
        /** It indicates whether the last push was a dummy push or not. */
        bool dummy;
        /** It indicates wheter the current states C and D are valid or not. 
          * There are valid if all of the ActiveTritAtCAndD of the unit-list 
          * have a neighbor by "constraint on D" and a neighbor by "constraint 
          * on C".
          */
        bool valid;
        /** It indicates wheter the current states C and D are a new valid pattern
          * or or not.
          */ 
        bool newValidPattern; 
    public:
        ActiveStatesCAndDCache(); 
        void pushDummy(); 
        void push(const ActiveTritAtCAndD& trit);
        void pop(const ActiveTritAtCAndD& trit);
        /** It uses the attribute nrTritsOnSameTryteColumnAtC to know if 
          * @trit has a neighbor by constraint at C in the unit-list. 
          */ 
        bool hasNeighborAtC(const ActiveTritAtCAndD& trit) const; 
        /** It uses the attribute nrTritsOnSameColumnAtD to know if 
          * @trit has a neighbor by constraint at D in the unit-list. 
          */ 
        bool hasNeighborAtD(const ActiveTritAtCAndD& trit) const;
    private:
        bool isAValidPattern(const ActiveTritAtCAndD& last) const;
};

/** This class is used to compute a lower bound on the weights of a 
  * 3-round trail core (A, B, C, D), where the activity patterns of
  * C and D are given by the active states stateC and stateD. For that 
  * it uses the attributes of the ActiveStatesCAndDCache.
  */
class KK_TrailCoreCostFunction {
    public: 
        KK_TrailCoreCostFunction(){}; 
        unsigned int getCost(const vector<ActiveTritAtCAndD>& unitList,
                             const ActiveStatesCAndDCache& cache) const;
}; 

/** This class is used to find all the 3-round trail cores (A, B, C, D)
  * with B and D in the kernel up to some weight. 
  */ 
class KK_TrailCores {
    private:
        /** The maximum weight allowed for the 3-round trail cores. */ 
        long double T3; 
        /** The name of the file that will contain the trail cores. */ 
        string file_KK_TrailCores; 
    public: 
        KK_TrailCores(long double T3); 
        void generate_KK_trailCores();
};

#endif 
