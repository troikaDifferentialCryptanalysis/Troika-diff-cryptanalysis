/* sbox.h */ 
#ifndef SBOX_H
#define SBOX_H 

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include "state.h"

/** Class used to store an input (resp. output) tryte difference compatible with 
  * an implicit and particular output (resp. input) difference of the S-box, 
  * denoted here diff.
  */
class TryteSTCompatible : public Tryte {
    public:
        /**  The weight of the transition. */
        Weight weight;
    public:
        TryteSTCompatible() : weight(0) {}
        /* The constructor. 
         * @param aTryte  Input (resp. output) tryte difference compatible 
         *                through the map SubTrytes with an implicit and
         *                particular output (resp.input) difference diff. 
         * @param aWeight Weight of the transition aValue --S--> diff 
         *                (resp. diff --S--> aValue ).
         */
        TryteSTCompatible(unsigned int aValue, Weight aWeight) 
        : Tryte(aValue), weight(aWeight){} 
        friend ostream & operator << (ostream& fout,
                                      const TryteSTCompatible& aTryte);

};

/** This class is used during an in-kernel backward extension.
  * (see Appendix E.1 Backward extension inside the kernel)
  * It stores a column of trytes' differences compatible through the map 
  * SubTrytes with an implicit and particular column of trytes'differences,
  * denoted here diff := |diff[0] | diff[1] | diff[2] |.
  */
class TryteColumnSTCompatible {
    public:
        /** For  0 ≤ y < 3, trytes[y] contains the tryte value of the tryte 
          * difference of coordinate y of the column of trytes.
          */
        Tryte trytes[3]; 
        /** The total Hamming weight of the column of trytes. */
        unsigned int hammingWeight; 
        /** The sum of the weights of the 3 transitions trytes[y] --S--> diff[y]
         *  (or diff[y] --S--> trytes[y])
         */
        Weight weight;
    public:
        /* The constructor. 
         * @param aTryte0 The tryte difference of coordinate y = 0. 
         *.@param aTryte1 The tryte difference of coordinate y = 1.
         * @param aTryte2 The tryte difference of coordinate y = 2.
         */
        TryteColumnSTCompatible(TryteSTCompatible aTryte0, 
                                TryteSTCompatible aTryte1, 
                                TryteSTCompatible aTryte2);
        /** An ordering operator. Elements are ordering according to their
          * their contribution to the total weight of the backward extension
          * (see Appendix E.1 Backward extension inside the kernel).
          * The cost of a column of trytes, denoted tryteCol, is defined by
          * c(tryteCol) := 2 * tryteCol.hammingWeight + tryteCol.weight
          * @param tryteColumn1 The first column of trytes to compare 
          * @param tryteColumn2 The second column of trytes to compare. 
          * @return true if c(tryteColumn1) < c(tryteColumn2), false otherwise.
          */
        friend bool operator < (const TryteColumnSTCompatible& tryteColumn1, 
                                const TryteColumnSTCompatible& tryteColumn2);
        friend ostream & operator << (ostream& fout,
             const TryteColumnSTCompatible& aTryteColumnSTCompatible);
};

#endif
