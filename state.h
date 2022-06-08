/* state.h */ 
#ifndef STATE_H
#define STATE_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <vector>
#include "types.h"

class TrytePosition;
class Tryte;
class TroikaState;

/** Class containing the x and z coordinates of a column. */
class ColumnPosition
{
    public:
        /** The x-coordinate of the column, 0 ≤ x < 9. */
        unsigned int x;
        /** The z-coordinate of the column, 0 ≤ z < 27. */
        unsigned int z;
    public:
        ColumnPosition() : x(0), z(0) {}
        ColumnPosition(unsigned int ax, unsigned int az) : x(ax), z(az) {}
        unsigned int getXplus9Z() const { return x + COLUMNS*z; }
        friend ostream & operator << (ostream& fout, const ColumnPosition& aCP);
};

/** Class containing the x, y and z coordinates of a trit. */
class TritPosition {
    public:
        /** The x-coordinate of the trit, 0 ≤ x < 9. */
        unsigned int x;
        /** The y-coordinate of the trit, 0 ≤ y < 3. */
        unsigned int y;
        /** The z-coordinate of the trit, 0 ≤ z < 27. */
        unsigned int z;
    public: 
        TritPosition() : x(0), y(0), z(0) {}
        TritPosition(unsigned int ax, unsigned int ay, unsigned int az);
        TritPosition(const TritPosition& t);
        void set(unsigned int ax=0, unsigned int ay=0, unsigned int az=0);
        /** It sets the coordinates of the first, second or third trit of
          * a tryte.  
          * @param aTrytePosition The position of the tryte. 
          * @param tritIndex      The index 0, 1 or 2 of the trit of the tryte. 
          */
        void set(TrytePosition aTrytePosition, unsigned int tritIndex); 
        /** It increments the coordinates (x, y) of the trit while maintaining 
          * the z coordinate. It first attempts to increment y. If it is not
          * possible, it resets y and increments x.
          * @return false if (x, y) = (8, 2), true otherwise.
          */
        bool setNextXY(); 
        /** It translates the position along the x axis.
          * @param  dx  The amount of translation (must be positive).
          */
        void xTranslate(unsigned int dx);
        /** It translates the position along the y axis.
          * @param  dy  The amount of translation (must be positive).
          */
        void yTranslate(unsigned int dy);
        /** It translates the position along the z axis.
          * @param  dz  The amount of translation (must be positive).
          */
        void zTranslate(unsigned int dz); 
        /** It applies the ShiftRows and ShiftLanes permutation (SRSL) to the coordinates of the trit. */
        void SRSL();
        /** It applies the inverse of ShiftLanes and the inverse of ShiftRows to the coordinates of the trit. */
        void invSRSL();
        /** @return The image of the trit's coordinates by the ShiftRows and 
          *         ShiftLanes permutations.
          */
        TritPosition getSRSL() const;
        /** @return The antecedent of the trit's coordinates by the ShiftRows
          *         and ShiftLanes permutations. 
          */
        TritPosition getInvSRSL() const; 
        unsigned int getXplus9Y() const { return x + COLUMNS*y; }
        bool operator == (const TritPosition& other) const;
        /** Arbitrary ordering operator. It is the lexicographic order [z, x, y].
         * @param  otherTrit The trit position at the right of the operator.
         * @return True iff this object comes before @other.
         */
        bool operator < (const TritPosition& other) const;
        friend ostream & operator << (ostream& fout, const TritPosition& t);
};

/** Class containing the x, y and z coordinates of a tryte. */
class TrytePosition
{
    public:
        /** The x-coordinate of the tryte, 0 ≤ x < 3. */
        unsigned int x; 
        /** The y-coordinate of the tryte, 0 ≤ y < 3. */
        unsigned int y; 
        /** The z-coordinate of the tryte, 0 ≤ z < 27. */
        unsigned int z;
    public: 
        TrytePosition() : x(0), y(0), z(0) {}
        TrytePosition(unsigned int ax, unsigned int ay, unsigned int az) : 
            x(ax), y(ay), z(az) {}  
        /** A constructor.  
          * @param t TritPosition used by the constructor to initialize the
          *          tryte with the coordinates of the tryte containing the
          *          trit @t.          
          */
        TrytePosition(const TritPosition& t) : x(t.x/3), y(t.y), z(t.z) {}
        bool operator == (const TrytePosition& other) const;
        /** Arbitrary ordering operator. It is the lexicographic order [z, x, y].
          * @param  otherTryte  The tryte position at the right of the operator.
          * @return True iff this object comes before @a otherTryte.
          */
        bool operator < (const TrytePosition& other) const;
        friend ostream & operator << (ostream &fout, const TrytePosition &trytePosition);
};

/** Class containing the value of a tryte. It is used to perform the 
  * coordinate-wise addition between two trytes to compute the differential
  * distribution table of the Troika S-box. 
  * The values of the trits and the hamming weight of the 27 trytes are stored 
  * in the array TRYTES. 
  */
class Tryte {
    public:
        /** The value of the tryte difference. It is an integer between 0 and 
          * 26.
          * | t0 | t1 | t2 | <--> 9 * t0 + 3 * t1 + t2. 
          * where ti is the trit of index i of the tryte.
          */
        unsigned int value; 
    public:
        Tryte() : value(0) {} 
        /** The constructor.
          * @param aValue Must be an integer between 0 and 26. 
          */ 
        Tryte(unsigned int aValue) : value(aValue) {}
        /** @param i  The index 0 ≤ i < 3 of the trit. 
          * @return   The value of the trit of index i; 
          */
        unsigned int trit(unsigned int i) const;
        /** @return The hamming weight of the tryte. */ 
        unsigned int getHammingWeight() const; 
        explicit operator unsigned int() const { return value; }
        /** The addition operator. It is done trit by trit. 
          * @param other The tryte at the right of the operator. 
          */ 
        Tryte operator + (const Tryte& other) const;
        /** The subtraction operator. It is done trit by trit. 
          * @param other The tryte at the right of the operator. 
          */ 
        Tryte operator - (const Tryte& other) const;
        bool operator == (const Tryte& other) const;
        friend ostream & operator << (ostream &fout, const Tryte &tryte);
}; 

/** Class used to represent the trit-activity pattern of a Troika state (it 
  * stores the positions of the active trits of the state).
  */
class ActiveState {
    public:
        /** Array used to store the positions of the active trits of the 
          * state, seen as 27 slices.
          * If the bit activeSlices[z] & (1 << x + 9*y) equals 1, that means 
          * that the trit of x-coordinate, y-coordinate and z-coordinate is
          * active.
          */
        UINT32 activeSlices[SLICES]; 
    public:
        ActiveState(); 
        /** A constructor. It initializes the active trits' positions of the
          * active state at the same positions as those of the TroikaState.
          * @param  state  The TroikaState used for the initialization.
          */ 
        ActiveState(const TroikaState& state); 
        /** It sets the active trits' positions at the same positions as those
          * of the TroikaState. 
          * @param state The TroikaState used to set the trit-activity pattern.
          */
        void set(const TroikaState& state);  
        void activateTrit(unsigned int x, unsigned int y, unsigned int z); 
        void activateTrit(const TritPosition& t); 
        void deactivateTrit(unsigned int x, unsigned int y, unsigned int z); 
        void deactivateTrit(const TritPosition &t);
        /** It activates or deactivates the trit of position @t
          * @param value True if the trit must be activated, false otherwise. 
          * @param t     The position of the trit
          */
        void setTritValue(bool value, const TritPosition& t); 
        bool isTritActive(unsigned int x, unsigned int y, unsigned int z) const;
        bool isTritActive(const TritPosition& t) const;
        bool isTheTritInAnActiveTryte(unsigned int x, unsigned int y, unsigned int z) const;
        bool isTheTritInAnActiveTryte(const TritPosition& t) const;
        bool isTryteActive(const TrytePosition& pos) const; 
        bool isTryteActive(int xTryte, int y, int z) const; 
        bool isSliceActive(unsigned int z) const; 
        bool isColumnActive(unsigned int x, unsigned int z) const; 
        unsigned int getNrActiveTrytes() const;
        /** It gets the active trits' positions of a column.
          * @param x The x-coordinate, 0 ≤ x < 9, of the column  
          * @param z The z-coordinate, 0 ≤ z < 27, of the column.
          * @return  An integer isYiActive such that, for  0 ≤ y < 3,
          *          (isYiActive & (1 << y) ) != 0 iff the trit of the column
          *          of coordinate y is active. 
          */
        UINT8 getIsYiActive(unsigned int x, unsigned int z) const; 
        /** It gets the value of the active tryte. It is the integer 
          * activeTryteValue such that  for  0 ≤ i < 3,
          * (activeTryteValue & (1 << i)) != 0 iff the trit of 
          * index i is active. 
          * @param xTryte The x-coordinate, 0 ≤ x < 3, of the active tryte.   
          * @param y The y-coordinate, 0 ≤ y < 3, of the active tryte. 
          * @param z The z-coordinate, 0 ≤ z < 27, of the active tryte. 
          * @return  The value of the active trit. 
          *          | t0 | t1 | t2 | <--> t0 + 2 * t1 + 4*t2
          *           where ti is the active trit of index i.
          */
        int getActiveTryte(unsigned int xTryte, unsigned int y, unsigned int z) const;
        /** It sets in the variable @biggestRepresentative the biggest
          * representative among the translated version of the state along the
          * z axis. The order relation compares the slices one after the other,
          * using a lexicographic order.
          * @param biggestRepresentative The state that is going to contain the
          *                              biggest representative. 
          */
        void setTheBiggestRepresentative(ActiveState& biggestRepresentative) const;
        bool operator == (const ActiveState& other) const;
        /** An arbitrary ordering operator. It compares the slices one after the other, 
          * using a lexicograohic order.
          * @param  otherState  The active state at the right of the operator.
          * @return True iff this object comes before @a otherState.
          */
        bool operator < (const ActiveState& other) const;
        void save(ostream& fout) const;      
        void load(istream& fin);
        friend ostream & operator << (ostream& fout, const ActiveState& aActiveState);
};

/* Class containing the values of the 27 trits of a lane. */
class TroikaLane {
    public:
        /** It indicates where the trits equal to 1 are. For 0 ≤ z < 27,
          * if the bit lane_1 & (1 << z) equals 1, that means that the trit
          * of z-coordinate z equals 1. 
          */
        UINT32 lane_1; 
        /** It indicates where the trits equal to 2 are. For 0 ≤ z < 27,
          * If the bit lane_2 & (1 << z) equals 1, that means that the trit
          * of z-coordinate z equals 2. 
          */
        UINT32 lane_2;
    public:
        TroikaLane(): lane_1(0), lane_2(0) {}
        void set_0(unsigned int z); 
        void set_1(unsigned int z); 
        void set_2(unsigned int z);
        void add_1(unsigned int z); 
        void add_2(unsigned int z);
        /** It translates the lane along the z axis.
          * @param  dx  The amount of translation (0 ≤ dz < 27).
          */
        void shiftBy(unsigned int z);
        friend TroikaLane operator+ (const TroikaLane& a, const TroikaLane& b);
        friend TroikaLane operator- (const TroikaLane& a, const TroikaLane& b); 
        friend ostream & operator << (ostream& fout, const TroikaLane& aLane);
};

/** Class containing the values of the 27 * 9 trits of a plane. The plane is 
  * seen as 9 lanes. It is used to store the parity plane and the theta-effect 
  * plane of the states during the generation of parity bare states. 
  */
class TroikaPlane {
    public: 
        TroikaLane lanes[COLUMNS];
    public: 
        TroikaPlane();
        void set_0(unsigned int x, unsigned int z) { lanes[x].set_0(z);} 
        void set_1(unsigned int x, unsigned int z) { lanes[x].set_1(z);}  
        void set_2(unsigned int x, unsigned int z) { lanes[x].set_2(z);}
        void addTritValue(unsigned int value, unsigned int x, unsigned int z);
        void multiplyTritBy2(unsigned int x, unsigned int z);
        bool isTritActive(unsigned int x, unsigned int z) const; 
        unsigned int getTrit(unsigned int x, unsigned int z) const;
        friend ostream & operator << (ostream &fout, const TroikaPlane &aPlane);
};

/* Class containing the 27 lanes of a state. */
class TroikaState {
    public:
        /* The 27 lanes of the state : 
         * lanes[x + 9 y] for  0 ≤ x < 9 and 0 ≤ y < 3.
         */
        TroikaLane lanes[COLUMNS * ROWS];
    public: 
        /* The default constructor. It initializes a state with no active trit. */
        TroikaState(){}
        /** It deactivates all the active trits of the state. */
        void setEmptyState();
        /** It applies the ShiftRows and ShiftLanes permutations to the state. */
        void SRSL(); 
        /** It applies the inverses of the ShiftLanes and the inverse of the ShiftRows permutations to the state. */
        void invSRSL();
        /** It applies the addColumnParity permutation to the state. */
        void addColumnParity(); 
        /** It applies the inverse of addColumnParity of the state. */
        void invAddColumnParity();
        /** It applies the ShiftRows, the ShiftLanes and the addColumnParity permutations to the state. */
        void L(); 
        /** It applies the inverses of the addColumnParity, ShiftLanes and ShiftRows permutations to the state. */
        void invL();
        /** It sets the image of @state by the inverses of the ShiftLanes and ShiftRows permutations. */
        void setInvSRSL(const TroikaState& state);
        /** It sets the image of @state by the inverses of the addColumnParity, ShiftLanes and ShiftRows permutations. */
        void setInvL(const TroikaState& state); 
        bool isInKernel() const; 
        int getTrit(unsigned int x, unsigned int y, unsigned int z) const;
        int getTrit(const TritPosition& t) const; 
        bool isTritActive(unsigned int x, unsigned int y, unsigned int z) const;
        bool isTritActive(const TritPosition& t) const;
        void set_0(const TritPosition& t);
        void setTritValue(unsigned int value, const TritPosition& t);
        void setTritValue(unsigned int value, unsigned int x, unsigned int y, unsigned int z);  
        void add_1(const TritPosition& t); 
        void add_2(const TritPosition& t);
        void addTritValue(unsigned int value, const TritPosition& t);
        void multiplyTritBy2(const TritPosition& t); 
        bool isTheTritInAnActiveTryte(const TritPosition& t) const; 
        bool isTryteActive(const TrytePosition& pos) const;
        bool isTryteActive(unsigned int xTryte, unsigned int y, unsigned int z) const;
        bool isSliceActive(unsigned int z) const; 
        unsigned int getNrActiveTrytes() const;
        unsigned int getNrActiveTrytesOfSlice(unsigned int z) const;
        /** It returns the minimum reverse weight of the minimum weight of the
          * state, that is twice the number of active trytes.
          */ 
        Weight getWeight() const;
        void save(ostream& fout) const;
        void load(istream& fin);
        /** It returns in an array the value of the column of coordinates
          * @x and @z. For 0 ≤ y < 3, the element of the array of index y 
          * corresponds to the trit of the state of coordinates (@x, y, @z). 
          */
        array<int, 3> getColumn(unsigned int x, unsigned int z) const; 
        Tryte getTryte(unsigned int xTryte, unsigned int y, unsigned int z) const;
        void setTryte(unsigned int xTryte, unsigned int y, unsigned int z, const Tryte& tryte);
        /** It translates the state along the z axis.
          * @param  dz  The amount of translation (must be positive).
          */
        void translate(unsigned int dz);
        bool operator == (const TroikaState& otherState) const; 
        bool operator != (const TroikaState& otherState) const; 
        /** An arbitrary ordering operator. The order relation is the lexicographic order
          * [lanes[0].lane_1, lanes[0].lane_2, ..., lanes[26].lane_1, lanes[26].lane_2].
          */
        bool operator < (const TroikaState& otherState) const;
        friend ostream & operator << (ostream &fout, const TroikaState &aState); 
};

#endif
