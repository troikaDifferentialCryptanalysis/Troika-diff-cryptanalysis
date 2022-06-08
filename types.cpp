/* types.cpp */
#include "types.h"

ostream & operator << (ostream& fout, const Weight& aWeight)
{
    if (aWeight.logPart) {
        fout << aWeight.integer << " + " << aWeight.logPart << " LOG "; 
        fout << " = " << aWeight.integer + LOG * aWeight.logPart; 
    } else {
        fout << aWeight.integer; 
    }
    return fout; 
}
