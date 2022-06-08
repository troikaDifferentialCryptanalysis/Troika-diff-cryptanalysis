/* sbox.cpp */ 
#include "sbox.h"

const UINT8 SBOX[27] = {6, 25, 17, 5, 15, 10,  4, 20, 24, 
                        0,  1,  2, 9, 22, 26, 18, 16, 14,
                        3, 13, 23, 7, 11, 12,  8, 21, 19};

ostream & operator << (ostream& fout, const TryteSTCompatible& aTryte)
{
    fout << (Tryte)aTryte << " : "; 
    fout << aTryte.weight; 
    return fout; 
}
