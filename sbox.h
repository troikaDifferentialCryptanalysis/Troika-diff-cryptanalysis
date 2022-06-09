/* sbox.h */ 
#ifndef SBOX_H
#define SBOX_H 

/*

This code uses functions and ideas from KeccakTools
(https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

The classes of this file deal with the differentials over the nonlinear
SubTrytes (abbreviated ST) map. They are used to see wether differentials are
valid or not and to know the weight of a valid differential over Subtrytes. 

*/

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
  * (see Appendix C.1 Backward extension inside the kernel)
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
          * (see Appendix C.1 Backward extension inside the kernel).
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

/** Class used to store : 
  * - the DDT of the S-box, 
  * - the input differences compatible through SubTrytes(ST) with a particular 
  *   output difference,
  * - the output differences compatible through SubTrytes (ST) with a particular
  *   input difference,
  * - the in-kernel columns of trytes that are compatible through with a
  *   particular output column of trytes (see Apendix C.1, backward extensions 
  *   inside the kernel). 
  */
class Sbox
{
    public: 
        /** The differential distribution table.
          * DDT[output][input] = # (x | S(x + input) - S(x) = output) ) */
        static unsigned int DDT[27][27];
        /** For  0 ≤ inputDiff < 27, outputDiff[inputDiff] is a vector 
          * containing all the output tryte differences (value of the difference 
          * and weight of the transition) compatible with inputDiff.
          */
        static vector<TryteSTCompatible> outputDiff[27];
        /** For  0 ≤ outputDiff < 27, inputDiff[outputDiff] is a vector 
          * containing all the input tryte differences (value of the difference 
          * and weight of the transition) compatible with outputDiff.
          */
        static vector<TryteSTCompatible> inputDiff[27];
        /** This attribute is used for backward in kernel extensions. 
          * For  0 ≤ tryteA2 ≤ tryteA1 ≤ tryteA0 < 27,
          * inKernelTryteColumnBeforeST[tryteA0][tryteA1][tryteA2]
          * contains a vector of columns of trytes of form 
          * | tryteB0 | tryteB1 | tryteB2 | such that the input difference
          * tryteBi is compatible by the Sbox with the output difference tryteAi
          * and the column of tryte | tryteB0 | tryteB1 | tryteB2 | is in the
          * kernel. The elements of that vector are sorted in ascending order
          * according to their cost (see the class TryteColumnSTCompatible).
          */
        static vector<vector<vector<vector<TryteColumnSTCompatible>>>>
                                    inKernelTryteColumnBeforeST;
    private: 
        static bool _init; 
    public: 
        Sbox(){}
        /** It checks wether a Troika input difference is compatible with a
          * Troika output difference by ST or not. If that is the case, it sets
          * the weight of the transition in @weightToSet. 
          * @param inputDifference  The Troika input difference. 
          * @param outputDifference The Troika output difference.
          * @param positionsForSTCompatibility Vector containing all the
          *                                    positions of the active trytes
          *                                    of inputDifference 
          *                                    /outputDifference.
          * @param weightToSet If @inputDifference is compatible with 
          *                    @outputDifference, the weight of the transition
          *                    is set in that variable. 
          * @return true if the input difference is compatible through
          *         ST with the output difference, false otherwise. 
          */
        bool areSTCompatible(const TroikaState& inputDifference, 
                const TroikaState& outputDifference, 
                const vector<TrytePosition>& positionsForSTCompatibility,
                Weight& weightToSet) const;
        /** Same as above, excepting that the position of the active trytes 
          * of inputDifference and outputDifference are not known in advance.
          */ 
        bool areSTCompatible(const TroikaState &inputDifference, 
                             const TroikaState &outputDifference,
                             Weight &weightToSet) const;
        void displayDDT() const; 
        void displayOutputCompatibleWith(const Tryte &input) const; 
        void displayInputCompatibleWith(const Tryte &output) const; 
    private:
        static bool init();
        static void DDTInitialization(); 
        static void outputInputDiffInitialization(); 
        static bool isInKernel(Tryte aTryte1, Tryte aTryte2, Tryte aTryte3);
        static vector<TryteColumnSTCompatible> getTryteColumnsBeforeST(
                                                                unsigned i, 
                                                                unsigned j, 
                                                                unsigned k);
        static void inKernelTryteColumnBeforeSTInitialization();
        bool areSTCompatible(Tryte inputTryte, Tryte outputTryte, 
                             Weight &accumulatedWeight) const;
}; 
#endif
