/* state.cpp */ 
#include "state.h"
#include "sbox.h"

static const int SHIFT_ROWS_PARAM[3] = {0, 1, 2};
                                    
static const int SHIFT_LANES_PARAM[27] =
/*     x=0   x=1   x=2   x=3   x=4   x=5   x=6   x=7   x=8       */
      {19   ,13   ,21   ,10   ,24   ,15   ,2    ,9    ,3    , // y=0
       14   ,0    ,6    ,5    ,1    ,25   ,22   ,23   ,20   , // y=1
       7    ,17   ,26   ,12   ,8    ,18   ,16   ,11   ,4  };  // y=2


ostream & operator << (ostream& fout, const ColumnPosition& aCP)
{
    fout << "(" << dec << aCP.x << ",-," << setw(2) << aCP.z << ") "; 
    return fout;
}

TritPosition::TritPosition(unsigned int ax, unsigned int ay, unsigned int az)
: x(ax), y(ay), z(az)
{ 
}

TritPosition::TritPosition(const TritPosition& t)
: x(t.x) , y(t.y) , z(t.z)
{ 
}

void TritPosition::set(unsigned int ax, unsigned int ay, unsigned int az)
{
    x = ax;
    y = ay;
    z = az;
}

void TritPosition::set(TrytePosition aTrytePosition, unsigned int tritIndex)
{
    x = 3 * aTrytePosition.x + tritIndex; 
    y = aTrytePosition.y; 
    z = aTrytePosition.z; 
}

bool TritPosition::setNextXY()
{
    if (y < 2) {
        y++;
        return true; 
    }
    if (x < 8) {
        y = 0; 
        x++; 
        return true; 
    }
    return false;
}

void TritPosition::xTranslate(unsigned int dx)
{
    x = (x + dx) % COLUMNS;  
}

void TritPosition::yTranslate(unsigned int dy)
{
    y = (y + dy) % ROWS; 
}

void TritPosition::zTranslate(unsigned int dz)
{
    z = (z + dz)  % SLICES; 
}

void TritPosition::SRSL()
{
    x = (x + 3 * SHIFT_ROWS_PARAM[y]) % COLUMNS; 
    z = (z + SHIFT_LANES_PARAM[getXplus9Y()]) % SLICES; 
}

void TritPosition::invSRSL()
{
    z = (z - SHIFT_LANES_PARAM[getXplus9Y()] + SLICES) % SLICES; 
    x = (x - 3 * SHIFT_ROWS_PARAM[y] + COLUMNS) % COLUMNS;
}

TritPosition TritPosition::getSRSL() const
{
    TritPosition t = *this; 
    t.SRSL(); 
    return t; 
}

TritPosition TritPosition::getInvSRSL() const
{
    TritPosition t = *this;
    t.invSRSL(); 
    return t; 
}

bool TritPosition::operator == (const TritPosition &other) const
{
    return (x == other.x && y == other.y && z == other.z);
}

bool TritPosition::operator < (const TritPosition &other) const
{
    if (z < other.z)
        return true; 
    else if (z == other.z) {
        if (x < other.x)
            return true;
        else if (x == other.x) {
            if (y < other.y) {
                return true; 
            }
        }
    }
    return false; 
}

ostream & operator << (ostream& fout, const TritPosition& t)
{
    fout << "(" << dec << t.x << "," << t.y << ","<< setw(2) << t.z << ") ";
    return fout;
}

bool TrytePosition::operator==(const TrytePosition &other) const
{
    return (x==other.x && y==other.y && z==other.z);
}

bool TrytePosition::operator<(const TrytePosition &other) const
{
    if (z < other.z)
        return true; 
    else if (z == other.z) {
        if (x < other.x)
            return true;
        else if (x == other.x) {
            if (y < other.y) {
                return true; 
            }
        }
    }
    return false; 
}

ostream & operator << (ostream& fout, const TrytePosition &trytePosition)
{
    fout << "(" << dec << trytePosition.x ; 
    fout << "," << trytePosition.y ; 
    fout << ","<< setw(2) << trytePosition.z << ") ";
    return fout;
}

/** For 0 ≤ i < 27 : 
  *  - TRYTES[i][0] contains the value of the trit of index 0 of the tryte i. 
  *  - TRYTES[i][1] contains the value of the trit of index 1 of the tryte i.  
  *  - TRYTES[i][2] contains the value of the trit of index 2 of the tryte i. 
  *  - TRYTES[i][3] contains the hamming weigt of the tryte. 
  */  
static const unsigned int TRYTES[27][4] = { {0, 0, 0, 0}, 
                                            {0, 0, 1, 1},
                                            {0, 0, 2, 1}, 
                                            {0, 1, 0, 1}, 
                                            {0, 1, 1, 2},
                                            {0, 1, 2, 2}, 
                                            {0, 2, 0, 1}, 
                                            {0, 2, 1, 2},
                                            {0, 2, 2, 2},
                                            {1, 0, 0, 1}, 
                                            {1, 0, 1, 2},
                                            {1, 0, 2, 2}, 
                                            {1, 1, 0, 2}, 
                                            {1, 1, 1, 3},
                                            {1, 1, 2, 3}, 
                                            {1, 2, 0, 2}, 
                                            {1, 2, 1, 3},
                                            {1, 2, 2, 3},
                                            {2, 0, 0, 1}, 
                                            {2, 0, 1, 2},
                                            {2, 0, 2, 2}, 
                                            {2, 1, 0, 2}, 
                                            {2, 1, 1, 3},
                                            {2, 1, 2, 3}, 
                                            {2, 2, 0, 2}, 
                                            {2, 2, 1, 3},
                                            {2, 2, 2, 3} };

unsigned int Tryte::trit(unsigned int i) const 
{
    return TRYTES[value][i];
}

unsigned int Tryte::getHammingWeight() const
{
    return TRYTES[value][3]; 
} 

Tryte Tryte::operator + (const Tryte& other) const
{
    unsigned int res = 0; 
    for (unsigned int i = 0; i < 3; i++) 
        res +=  pow(3, 2 - i) * (( trit(i) + other.trit(i) ) % 3);  ; 
    return Tryte(res); 
}

Tryte Tryte::operator - (const Tryte& other) const
{
    unsigned int res = 0; 
    for (unsigned int i = 0; i < 3; i++) 
        res +=  pow(3, 2 - i) * (( trit(i) - other.trit(i) + 3) % 3);  ; 
    return Tryte(res);
}

bool Tryte::operator == (const Tryte &otherTryte) const
{
    return (value == otherTryte.value);

}

ostream & operator << (ostream& fout, const Tryte& tryte)
{
    for (int i = 0; i < 3; i++)
        fout << tryte.trit(i); 
    return fout; 
}

ActiveState::ActiveState()
{ 
    for (unsigned z = 0; z < SLICES; z++)
        activeSlices[z] = 0; 
}

ActiveState::ActiveState(const TroikaState& state)
{ 
    for (unsigned z = 0; z < SLICES; z++)
        activeSlices[z] = 0; 

    for (unsigned int y = 0; y < ROWS; y++) {
        for (unsigned int x = 0; x < COLUMNS; x++) {
            unsigned int i = 9 * y + x; 
            if ( state.lanes[i].lane_1 | state.lanes[i].lane_2 ) { 
                for (unsigned int z = 0; z < SLICES; z++) {
                    if ( (1 << z ) & (state.lanes[i].lane_1 | state.lanes[i].lane_2 )) {
                        activateTrit(x, y, z); 
                    }   
                }
            }
        }
    }
}

void ActiveState::set(const TroikaState &state)
{
    *this = ActiveState(state);
}

inline UINT32 getSlicePoint(unsigned int x, unsigned int y)
{
    return (UINT32) 1 << (x + COLUMNS * y); 
}

void ActiveState::activateTrit(unsigned int x, unsigned int y, unsigned int z)
{
    activeSlices[z] = activeSlices[z] | getSlicePoint(x, y); 
}

void ActiveState::activateTrit(const TritPosition& t)
{
    activeSlices[t.z] = activeSlices[t.z] | getSlicePoint(t.x, t.y); 
}

void ActiveState::deactivateTrit(unsigned int x, unsigned int y, unsigned int z)
{
    activeSlices[z] = activeSlices[z] & (~getSlicePoint(x, y)); 
}

void ActiveState::deactivateTrit(const TritPosition &t)
{
    activeSlices[t.z] = activeSlices[t.z] & (~getSlicePoint(t.x, t.y)); 
}

void ActiveState::setTritValue(bool value, const TritPosition& t)
{
    if (value)
        activeSlices[t.z] = activeSlices[t.z] | getSlicePoint(t.x, t.y); 
    else 
        activeSlices[t.z] = activeSlices[t.z] & (~getSlicePoint(t.x, t.y)); 
}

bool ActiveState::isTritActive(unsigned int x, unsigned int y, unsigned int z) const
{
    return bool(activeSlices[z] & getSlicePoint(x, y)); 
}

bool ActiveState::isTritActive(const TritPosition& t) const
{
    return bool(activeSlices[t.z] & getSlicePoint(t.x, t.y)); 
}

/** @param  index  The 0 ≤ index < 27 of the trit of the slice.
  * @return The 3 bits of the tryte containing the trit of index @index. 
  */
inline UINT32 getSliceTryteFromTritIndex(unsigned int index)
{
    return 0x7 << (3*( index / 3));  
}

bool ActiveState::isTheTritInAnActiveTryte(unsigned int x, unsigned int y, unsigned int z) const
{
    return bool(activeSlices[z] & getSliceTryteFromTritIndex(x + 9 *y));
}

bool ActiveState::isTheTritInAnActiveTryte(const TritPosition& t) const
{
    return bool(activeSlices[t.z] & getSliceTryteFromTritIndex(t.getXplus9Y()));
}

bool ActiveState::isTryteActive(const TrytePosition &pos) const
{ 
    return bool ( 0x7 & (activeSlices[pos.z] >>  (3*pos.x + 9 * pos.y)) ); 
}

bool ActiveState::isTryteActive(int xTryte, int y, int z) const
{
    return bool ( 0x7 & (activeSlices[z] >>  (3*xTryte + 9 * y)) ); 
}

bool ActiveState::isSliceActive(unsigned int z) const 
{
    return bool (activeSlices[z]); 
}

bool ActiveState::isColumnActive(unsigned int x, unsigned int z) const
{
    for (unsigned int y = 0; y < ROWS; y++) {
        if (isTritActive(x, y, z))
            return true; 
    }
    return false; 
}

unsigned int ActiveState::getNrActiveTrytes() const
{
    unsigned int nrActiveTrytes = 0; 
    for (unsigned int z = 0; z < SLICES; z++) {
        for (unsigned int xTryte = 0; xTryte < 3; xTryte++) {
            for (unsigned int y = 0; y < ROWS; y++) {
                if (isTryteActive(xTryte, y, z))
                    nrActiveTrytes += 1; 
            }
        }
    }
    return nrActiveTrytes; 
}

UINT8 ActiveState::getIsYiActive(unsigned int x, unsigned int z) const
{
    UINT8 isYiActive = 0; 
    for (unsigned int y = 0; y < ROWS; y++) {
        if (isTritActive(x, y, z))
            isYiActive ^= (1 << y); 
    } 
    return isYiActive; 
}

int ActiveState::getActiveTryte(unsigned int xTryte, 
                                         unsigned int y, 
                                         unsigned int z) const
{
    return 0x7 & (activeSlices[z] >>  (3*xTryte + 9 * y)); 
}

void ActiveState::setTheBiggestRepresentative(ActiveState& biggestRepresentative) const
{
    unsigned int bestDz = 0;
    unsigned int z; 

    for (unsigned int dz = 1; dz < SLICES; dz++) {
        z = 0; 
        while ( (z < SLICES - 1) 
                && (activeSlices[(z + bestDz) % SLICES] == activeSlices[(z + dz) % SLICES]))
            z++;
        if (activeSlices[(z + bestDz) % SLICES] < activeSlices[(z + dz) % SLICES])
            bestDz = dz;
    }
    for (unsigned int z = 0; z < SLICES; z++) {
        biggestRepresentative.activeSlices[z] = activeSlices[(z + bestDz) % SLICES]; 
    }
}

bool ActiveState::operator == (const ActiveState& other) const
{
    for (unsigned z = 0; z < SLICES; z++) {
        if (activeSlices[z] != other.activeSlices[z])
            return false; 
    }
    return true; 
}

bool ActiveState::operator < (const ActiveState& otherState) const
{
    for (unsigned z = 0; z < SLICES; z++) {
        if (activeSlices[z] > otherState.activeSlices[z])
            return false;
        if (activeSlices[z] < otherState.activeSlices[z])
            return true; 
    }
    return false; 
}

void ActiveState::save(ostream &fout) const
{
    fout << hex; 
    for (unsigned int z = 0; z < SLICES; z++) {
        fout << activeSlices[z] << " "; 
    }
    fout << endl; 
}

void ActiveState::load(istream &fin)
{
    fin >> hex; 
    for (unsigned int z = 0; z < SLICES; z++) {
        fin >> activeSlices[z];
    } 
}

ostream & operator << (ostream &fout, const ActiveState &aActiveState)
{
    for (unsigned int z = 0; z < SLICES; z++) {
        if (aActiveState.activeSlices[z] != 0) {
            fout <<"#### Slice " << z << "#### \n";
            for (unsigned int index = 0; index < COLUMNS*ROWS; index++) {
                if (aActiveState.activeSlices[z] & (1 << index))
                    fout << "1";
                else 
                    fout << "0"; 
                if (index % COLUMNS == 8)
                    fout << "\n"; 
                else 
                    fout << " "; 
            }
        }
    }
    return fout; 
}

inline UINT32 _1(TroikaLane lane)
{
    return lane.lane_1; 
}

inline UINT32 _2(TroikaLane lane)
{
    return lane.lane_2; 
}

inline UINT32 _0(TroikaLane lane)
{
    return (~lane.lane_1) & (~lane.lane_2); 
}

void TroikaLane::set_0(unsigned int z)
{
    lane_1 = lane_1 & (~(1 << z));
    lane_2 = lane_2 & (~(1 << z)); 
}

void TroikaLane::set_1(unsigned int z)
{
    lane_1 = lane_1 | (1 << z);
    lane_2 = lane_2 & (~(1 << z)); 
}

void TroikaLane::set_2(unsigned int z)
{
    lane_1 = lane_1 & (~(1 << z));
    lane_2 = lane_2 | (1 << z); 
}

void TroikaLane::add_1(unsigned int z)
{
    UINT32 temp = lane_1; 
    lane_1 = (lane_1 & ~(1 << z)) | ( ~lane_1 & ~lane_2 & (1 << z) ); 
    lane_2 = (lane_2 & ~(1 << z)) | (temp & (1 << z)); 
}

void TroikaLane::add_2(unsigned int z)
{
    UINT32 temp = lane_2; 
    lane_2 = ( lane_2 & ~(1 << z) ) | (~lane_1 & ~lane_2 & (1 << z)); 
    lane_1 = ( lane_1 & ~(1 << z) ) | (temp & (1 << z));  
}

void TroikaLane::shiftBy(unsigned int z)
{
    lane_1 = ((lane_1 << z) | (lane_1 >> (27-z))) & 0x07ffffff;
    lane_2 = ((lane_2 << z) | (lane_2 >> (27-z))) & 0x07ffffff;
}

TroikaLane operator+ (const TroikaLane &a, const TroikaLane &b)
{
    TroikaLane result;
    result.lane_1 = (_1(a) & _0(b)) | (_0(a) & _1(b)) | (_2(a) & _2(b)); 
    result.lane_2 = (_2(a) & _0(b)) | (_0(a) & _2(b)) | (_1(a) & _1(b)); 
    return result; 
}

TroikaLane operator- (const TroikaLane &a, const TroikaLane &b)
{
    TroikaLane result;
    result.lane_1 = (_1(a) & _0(b)) | (_2(a) & _1(b)) | (_0(a) & _2(b)); 
    result.lane_2 = (_2(a) & _0(b)) | (_1(a) & _2(b)) | (_0(a) & _1(b)); 
    return result; 
}

ostream & operator << (ostream& fout, const TroikaLane& aLane)
{
    for (unsigned int z = 0; z < 27; z++) {
        if (aLane.lane_1 & (1 << z))
            fout << "1 "; 
        else if (aLane.lane_2 & (1 << z))
            fout << "2 "; 
        else 
            fout << "0 "; 
    }
    return fout; 
}

TroikaPlane::TroikaPlane()
{
    for (unsigned int x = 0; x < COLUMNS; x++) {
        lanes[x].lane_1 = 0;
        lanes[x].lane_2 = 0; 
    }
}

void TroikaPlane::addTritValue(unsigned int value, unsigned int x, unsigned int z)
{
    value = value % 3; 
    TroikaLane toAdd; 
    if (value % 3 == 1)
        toAdd.lane_1 = (1 << z); 
    if (value % 3 == 2)
        toAdd.lane_2 = (1 << z); 
    lanes[x] = lanes[x] + toAdd; 
}    

void TroikaPlane::multiplyTritBy2(unsigned int x, unsigned int z) 
{
    // FAUX !
    UINT32 mask = (lanes[x].lane_1 ^ lanes[x].lane_2) & (1 << z); 
    lanes[x].lane_1 ^= mask;
    lanes[x].lane_2 ^= mask;  
}

bool TroikaPlane::isTritActive(unsigned int x, unsigned int z) const
{
    return bool (  (lanes[x].lane_1 & (1 << z))
                 | (lanes[x].lane_2 & (1 << z)) );  
}

unsigned int TroikaPlane::getTrit(unsigned int x, unsigned int z) const
{
    if (lanes[x].lane_1 & (1 << z))
        return 1; 
    if (lanes[x].lane_2 & (1 << z))
        return 2;
    return 0; 
}

ostream & operator << (ostream &fout, const TroikaPlane &aPlane)
{
    UINT32 mask = (1 << 26); 
    for (unsigned int z = 0; z < SLICES; z++) {
        for (unsigned x = 0; x < COLUMNS; x++) {
            if (aPlane.lanes[x].lane_1 & mask) 
                fout << "1"; 
            else if (aPlane.lanes[x].lane_2 & mask) 
                fout << "2"; 
            else 
                fout << "0"; 
            fout << " ";     
        }
        mask = mask >> 1;
        fout << "\n"; 
    }
    return fout; 
}

void TroikaState::setEmptyState()
{
    for (unsigned int i = 0; i < 27; i++) {
        lanes[i].lane_1 = 0; 
        lanes[i].lane_2 = 0; 
    }
}

void TroikaState::SRSL()
{
    // shiftRows
    TroikaState temp_state;
    const int shifts[27] = { 0,  1,  2,  3,  4,  5,  6, 7, 8,
                            12, 13, 14, 15, 16, 17, 9, 10, 11,
                            24, 25, 26, 18, 19, 20,21, 22, 23};
    for (int i = 9; i < 27; i++) {
        temp_state.lanes[i] = lanes[i]; 
    }
    for (int i = 9; i < 27; ++i) {
        lanes[shifts[i]] = temp_state.lanes[i];
    }
    
    // shiftLanes
    for (unsigned i = 0; i < 27; i++) {
        lanes[i].shiftBy(SHIFT_LANES_PARAM[i]);
    }
}

void TroikaState::invSRSL()
{
    TroikaState temp_state;
    const int shifts[27] = { 0,  1,  2,  3,  4,  5,  6,  7, 8,
                            15, 16, 17,  9, 10, 11, 12, 13, 14, 
                            21, 22, 23, 24, 25, 26, 18, 19, 20};
    // invShiftLanes
    for (unsigned i = 0; i < 27; i++) {
        lanes[i].shiftBy(27 - SHIFT_LANES_PARAM[i]);
    }
    // invShiftRows
    for (int i = 9; i < 27; i++) {
        temp_state.lanes[i] = lanes[i]; 
    }
    for (int i = 9; i < 27; ++i) {
        lanes[shifts[i]] = temp_state.lanes[i];
    }
}

void TroikaState::addColumnParity() 
{
    TroikaLane parity[COLUMNS]; 
    TroikaLane sum, l1, l2; 

    // Compute parity
    for (unsigned int x = 0; x < COLUMNS; x++) {
        sum.lane_1 = 0; 
        sum.lane_2 = 0; 
        for (unsigned int y = 0; y < ROWS; y++) {
            sum = sum + lanes[x + COLUMNS * y]; 
        }
        parity[x] = sum; 
    }

    // Add parity
    for (unsigned int y = 0; y < ROWS; y++) {
        for (unsigned int x = 0; x < COLUMNS; x++) {
            l1 = parity[(x - 1 + COLUMNS) % COLUMNS];  
            l2 = parity[(x + 1) % COLUMNS]; 
            l2.shiftBy(SLICES - 1); 
            lanes[x + COLUMNS * y] = lanes[x + COLUMNS * y] + l1 + l2; 
        }
    }
}

void TroikaState::invAddColumnParity() 
{
    TroikaLane parity[COLUMNS]; 
    TroikaLane sum, l1, l2; 

    // Compute parity
    for (unsigned int x = 0; x < COLUMNS; x++) {
        sum.lane_1 = 0; 
        sum.lane_2 = 0; 
        for (unsigned int y = 0; y < ROWS; y++) {
            sum = sum + lanes[x + COLUMNS * y]; 
        }
        parity[x] = sum; 
    }

    // Subtract parity
    for (unsigned int y = 0; y < ROWS; y++) {
        for (unsigned int x = 0; x < COLUMNS; x++) {
            l1 = parity[(x - 1 + COLUMNS) % COLUMNS];  
            l2 = parity[(x + 1) % COLUMNS]; 
            l2.shiftBy(SLICES - 1); 
            lanes[x + COLUMNS * y] = lanes[x + COLUMNS * y] - l1 - l2; 
        }
    }
}

void TroikaState::L()
{
    SRSL(); 
    addColumnParity();
}

void TroikaState::invL()
{
    invAddColumnParity();
    invSRSL(); 
}

void TroikaState::setInvSRSL(const TroikaState& state)
{   
    *this = state;  
    invSRSL(); 
}

void TroikaState::setInvL(const TroikaState& state)
{   
    *this = state;
    invL();
}

bool TroikaState::isInKernel() const
{
    TroikaLane sum; 
    // Compute parity
    for (unsigned int x = 0; x < COLUMNS; x++) {
        sum.lane_1 = 0; 
        sum.lane_2 = 0; 
        for (unsigned int y = 0; y < ROWS; y++) {
            sum = sum + lanes[x + COLUMNS * y]; 
        }
        if (sum.lane_1 || sum.lane_2 )
            return false; 
    }
    return true; 
}

int TroikaState::getTrit(unsigned int x, unsigned int y, unsigned int z) const
{
    if (lanes[9 * y + x].lane_1 & (1 << z))
        return 1; 
    if (lanes[9 * y + x].lane_2 & (1 << z))
        return 2; 
    return 0; 
}
     
int TroikaState::getTrit(const TritPosition& t) const
{
    if (lanes[9 * t.y + t.x].lane_1 & (1 << t.z))
        return 1; 
    if (lanes[9 * t.y + t.x].lane_2 & (1 << t.z))
        return 2; 
    return 0; 
}

bool TroikaState::isTritActive(unsigned int x, unsigned int y, unsigned int z) const
{
    return bool   (   (lanes[9 * y + x].lane_1 & (1 << z)) 
                   |  (lanes[9 * y + x].lane_2 & (1 << z)));
}

bool TroikaState::isTritActive(const TritPosition& t) const
{
    return bool (  (lanes[9 * t.y + t.x].lane_1 & (1 << t.z)) 
                 | (lanes[9 * t.y + t.x].lane_2 & (1 << t.z)));
}


void TroikaState::set_0(const TritPosition& t)
{
    lanes[9 * t.y + t.x].set_0(t.z); 
}

void TroikaState::setTritValue(unsigned int value, const TritPosition& t)
{
    if (value == 0)
        lanes[9 * t.y + t.x].set_0(t.z);
    if (value == 1)
        lanes[9 * t.y + t.x].set_1(t.z);
    if (value == 2)
        lanes[9 * t.y + t.x].set_2(t.z); 
}

void TroikaState::setTritValue(unsigned int value, unsigned int x, unsigned int y, unsigned int z)
{
    if (value == 0)
        lanes[9 * y + x].set_0(z);
    if (value == 1)
        lanes[9 * y + x].set_1(z);
    if (value == 2)
        lanes[9 * y + x].set_2(z); 
}

void TroikaState::add_1(const TritPosition& t)
{
    lanes[9 * t.y + t.x].add_1(t.z); 
}

void TroikaState::add_2(const TritPosition& t)
{
    lanes[9 * t.y + t.x].add_2(t.z); 
}

void TroikaState::addTritValue(unsigned int value, const TritPosition& t)
{
    value = value % 3; 
    if (value == 0)
        return; 
    if (value == 1)
        lanes[9 * t.y + t.x].add_1(t.z); 
    else
        lanes[9 * t.y + t.x].add_2(t.z); 
}

void TroikaState::multiplyTritBy2(const TritPosition& t)
{
    UINT32 mask = (lanes[9 * t.y + t.x].lane_1 ^ lanes[9 * t.y + t.x].lane_2)
                  & (1 << t.z); 
    lanes[9 * t.y + t.x].lane_1 ^= mask; 
    lanes[9 * t.y + t.x].lane_2 ^= mask; 
}

bool TroikaState::isTheTritInAnActiveTryte(const TritPosition& t) const
{
    unsigned x = 3 * (t.x / 3);
    for (unsigned xOffset = 0; xOffset < 3; xOffset++) {
        if (getTrit(x + xOffset, t.y, t.z))
            return true; 
    }
    return false; 
}

bool TroikaState::isTryteActive(const TrytePosition& pos) const
{
    return bool   (   (lanes[9 * pos.y + 3 * pos.x].lane_1 & (1 << pos.z)) 
                   |  (lanes[9 * pos.y + 3 * pos.x].lane_2 & (1 << pos.z))
                   |  (lanes[9 * pos.y + 3 * pos.x + 1].lane_1 & (1 << pos.z)) 
                   |  (lanes[9 * pos.y + 3 * pos.x + 1].lane_2 & (1 << pos.z))
                   |  (lanes[9 * pos.y + 3 * pos.x + 2].lane_1 & (1 << pos.z)) 
                   |  (lanes[9 * pos.y + 3 * pos.x + 2].lane_2 & (1 << pos.z))
                   );
}


bool TroikaState::isTryteActive(unsigned int xTryte, unsigned int y, unsigned int z) const
{
    return bool   (   (lanes[9 * y + 3 * xTryte].lane_1 & (1 << z)) 
                   |  (lanes[9 * y + 3 * xTryte].lane_2 & (1 << z))
                   |  (lanes[9 * y + 3 * xTryte + 1].lane_1 & (1 << z)) 
                   |  (lanes[9 * y + 3 * xTryte + 1].lane_2 & (1 << z))
                   |  (lanes[9 * y + 3 * xTryte + 2].lane_1 & (1 << z)) 
                   |  (lanes[9 * y + 3 * xTryte + 2].lane_2 & (1 << z))
                   );

}

bool TroikaState::isSliceActive(unsigned int z) const
{
    for (unsigned int i = 0; i < 27; i++) {
        if (lanes[i].lane_1 & (1 << z))
            return true; 
        if (lanes[i].lane_2 & (1 << z))
            return true; 
    }
    return false; 
}

unsigned int TroikaState::getNrActiveTrytes() const
{
    unsigned int nrActiveTrytes = 0; 
    for (unsigned int z = 0; z < SLICES; z++) {
        for (unsigned int xTryte = 0; xTryte < 3; xTryte++) {
            for (unsigned int y = 0; y < ROWS; y++) {
                if (isTryteActive(xTryte, y, z))
                    nrActiveTrytes += 1; 
            }
        }
    }
    return nrActiveTrytes; 
}

unsigned int TroikaState::getNrActiveTrytesOfSlice(unsigned int z) const
{
    unsigned int nrActiveTrytes = 0; 
    for (unsigned int xTryte = 0; xTryte < 3; xTryte++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            if (isTryteActive(xTryte, y, z))
                nrActiveTrytes += 1; 
        }
    }
    return nrActiveTrytes;
}

Weight TroikaState::getWeight() const
{
    return Weight(2 * getNrActiveTrytes()); 
}

void TroikaState::save(ostream &fout) const
{
    fout << hex; 
    for (unsigned int i = 0; i < 27; i++) {
        fout << lanes[i].lane_1 << " "; 
        fout << lanes[i].lane_2 << " "; 
    }
}

void TroikaState::load(istream &fin)
{
    fin >> hex; 
    for (unsigned int i = 0; i < 27; i++) {
        fin >> lanes[i].lane_1; 
        fin >> lanes[i].lane_2; 
    }
}

array<int, 3> TroikaState::getColumn(unsigned int x, unsigned int z) const 
{
    array<int, 3> col; 
    for (unsigned int y = 0; y < 3; y++) {
        col[y] = getTrit(x, y, z); 
    }
    return col; 
}

Tryte TroikaState::getTryte(unsigned int xTryte, unsigned int y, unsigned int z) const
{
    unsigned int value = 0; 
    for (unsigned int tritIndex = 0; tritIndex < 3; tritIndex++) {
        value += getTrit(3 * xTryte + tritIndex, y, z) * pow(3, 2 - tritIndex); 
    }
    return Tryte(value); 
}

void TroikaState::setTryte(unsigned int xTryte, unsigned int y,
                           unsigned int z, const Tryte& tryte)
{
    setTritValue(tryte.trit(0), 3 * xTryte    , y, z); 
    setTritValue(tryte.trit(1), 3 * xTryte + 1, y, z); 
    setTritValue(tryte.trit(2), 3 * xTryte + 2, y, z); 
}

void TroikaState::translate(unsigned dz) 
{
    for (unsigned int i = 0; i < 27; i++) {
        lanes[i].shiftBy(dz);
    }
}

bool TroikaState::operator == (const TroikaState& otherState) const
{
    for (unsigned int i = 0; i < 27; i++) {
        if (lanes[i].lane_1 != otherState.lanes[i].lane_1)
            return false; 
        if (lanes[i].lane_2 != otherState.lanes[i].lane_2)
            return false; 
    }
    return true; 
}

bool TroikaState::operator != (const TroikaState& otherState) const
{
    return ! ( *this == otherState );
}

bool TroikaState::operator < (const TroikaState& otherState) const
{
    for (unsigned int i = 0; i < 27; i++) {
        if (lanes[i].lane_1 > otherState.lanes[i].lane_1)
            return false; 
        if (lanes[i].lane_1 < otherState.lanes[i].lane_1)
            return true; 
        if (lanes[i].lane_2 > otherState.lanes[i].lane_2)
            return false; 
        if (lanes[i].lane_2 < otherState.lanes[i].lane_2)
            return true; 
    }
    return false; 
}

ostream & operator << (ostream &fout, const TroikaState &aState)
{ 
    for (unsigned z = 0; z < SLICES; z++) { 
        bool isSliceActive = false;
        for (unsigned int y = 0; y < ROWS; y++) {
            if (isSliceActive)
                break; 
            for (unsigned int x = 0; x < COLUMNS; x++) {
                if (aState.getTrit(x, y, z)) {
                    isSliceActive = true; 
                    break; 
                }
            }
        }
        if (isSliceActive) {
            fout <<"\n#### Slice " << z << "#### \n";
            for (unsigned int y = 0; y < ROWS; y++) {
                for (unsigned int x = 0; x < COLUMNS; x++) {
                    fout << aState.getTrit(x, y, z) << " " ; 
                }
                fout << endl; 
            }
        }
    }
    return fout; 
}
