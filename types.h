/* types.h */
#ifndef TYPES_H
#define TYPES_H 
/*

This code uses functions and ideas from KeccakTools
(https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

*/

#include <iostream>
#include <assert.h>
#include <sstream>

using namespace std;
typedef unsigned char       UINT8;
typedef unsigned short      UINT16;
typedef unsigned int        UINT32;
typedef unsigned long long  UINT64;

#define COLUMNS  9 
#define ROWS     3
#define SLICES   27
#define DIAGONAL 9
#define LOG      2.369070246428542692 /* - log <sub> 3 <\sub> (2 / 27) */

/** Exception with a string expressing the reason.
  */
class Exception {
    public:
        string reason;
    public:
        Exception(): reason() {}
        Exception(const string& aReason): reason(aReason) {}
};

/** This class is used to store the weight of a trail core. Since an active 
  * sbox adds to the total weight 2, LOG or 3, such a weight is of the form 
  * integer + LOG * logPart, where integer and logPart are integers.
  */
class Weight {
    public:
        unsigned int integer;
        unsigned int logPart;
    public: 
        Weight(): integer(0), logPart(0) {} 

        Weight(unsigned int aInteger, unsigned int aLogPart = 0):
            integer(aInteger), logPart(aLogPart){}
 
        explicit operator double long() const
        {
            return integer + LOG * logPart;
        }

        Weight operator + (const Weight& other) const
        {
            return Weight(integer + other.integer, 
                          logPart + other.logPart);
        }

        Weight operator - (const Weight& other) const
        {
            assert(integer >= other.integer);
            assert(logPart >= other.logPart); 
            return Weight(integer - other.integer,
                          logPart - other.logPart);
        }

        void operator += (const Weight& other)
        {
            integer += other.integer; 
            logPart += other.logPart;
        }

        void operator -= (const Weight& other)
        {
            assert(integer >= other.integer);
            assert(logPart >= other.logPart); 
            integer -= other.integer; 
            logPart -= other.logPart;
        }

        Weight operator * (unsigned int scalar) const
        {
            return Weight(scalar * integer, scalar * logPart);
        }

        bool operator < (const Weight& other) const
        {
            return ((long double)*this < (long double)other);
        }

        bool operator != (const Weight& other) const
        {
            return (integer != other.integer || logPart != other.logPart);
        }

        friend ostream & operator << (ostream &fout, const Weight &aWeight); 
};

#endif 
