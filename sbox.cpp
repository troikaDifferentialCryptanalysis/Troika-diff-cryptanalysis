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

TryteColumnSTCompatible::TryteColumnSTCompatible(
                                 TryteSTCompatible aTryte0, 
                                 TryteSTCompatible aTryte1, 
                                 TryteSTCompatible aTryte2)
{
    trytes[0] = aTryte0.value; 
    trytes[1] = aTryte1.value; 
    trytes[2] = aTryte2.value;
    hammingWeight =   aTryte0.getHammingWeight()
                    + aTryte1.getHammingWeight()
                    + aTryte2.getHammingWeight();
    weight =   aTryte0.weight
             + aTryte1.weight
             + aTryte2.weight; 
}

bool operator < (const TryteColumnSTCompatible &tryteColumn1, 
                 const TryteColumnSTCompatible &tryteColumn2)
{
    if (Weight(2 * tryteColumn1.hammingWeight) + tryteColumn1.weight
        < Weight(2 * tryteColumn2.hammingWeight) + tryteColumn2.weight)
        return true; 
    return false; 
}

ostream & operator << (ostream &fout, 
        const TryteColumnSTCompatible &aTryteColumnSTCompatible)
{
    fout << "cost : " << 2* aTryteColumnSTCompatible.hammingWeight + (long double)aTryteColumnSTCompatible.weight << endl;
    fout << "weight: " << aTryteColumnSTCompatible.weight << endl; 
    for (int i = 0; i < 3; i++)
        fout << aTryteColumnSTCompatible.trytes[i]  << endl;
    return fout; 
}

unsigned int Sbox:: DDT[27][27];
vector<TryteSTCompatible> Sbox::outputDiff[27];
vector<TryteSTCompatible> Sbox::inputDiff[27];
vector<vector<vector<vector<TryteColumnSTCompatible>>>>
                                    Sbox::inKernelTryteColumnBeforeST;
bool Sbox::_init = Sbox::init();

bool Sbox::init()
{
    Sbox::DDTInitialization(); 
    Sbox::outputInputDiffInitialization(); 
    Sbox::inKernelTryteColumnBeforeSTInitialization(); 
    return true; 
}

bool Sbox::areSTCompatible(const TroikaState& inputDifference, 
                           const TroikaState& outputDifference, 
                           const vector<TrytePosition>& positionsForSTCompatibility,
                           Weight& weightToSet) const
{
    vector<TrytePosition>::const_iterator it;
    Tryte inputTryte;
    Tryte outputTryte; 

    weightToSet = Weight(0); 

    for (it = positionsForSTCompatibility.begin();
         it != positionsForSTCompatibility.end(); it++) {
        inputTryte  = inputDifference.getTryte(it->x, it->y, it->z); 
        outputTryte = outputDifference.getTryte(it->x, it->y, it->z); 
        if (!areSTCompatible(inputTryte, outputTryte, weightToSet)) {
            return false; 
        }
    }
    return true; 
}


bool Sbox::areSTCompatible(const TroikaState &inputDifference, 
                             const TroikaState &outputDifference,
                             Weight &weightToSet) const
{
    vector<TrytePosition> positionsForSTCompatibility; 
    for (unsigned int xTryte = 0; xTryte < 3; xTryte++) {
        for (unsigned int y = 0; y < ROWS; y++) {
            for (unsigned int z = 0; z < SLICES; z++) {
                if (inputDifference.isTryteActive(xTryte, y, z)) {
                    positionsForSTCompatibility.push_back(TrytePosition(xTryte, y, z));
                }
            }
        }
    }
    return areSTCompatible(inputDifference, outputDifference,
                           positionsForSTCompatibility, weightToSet);
}


void Sbox::displayDDT() const
{ 
    for (unsigned int output = 0; output < 27; output++) {
        for (unsigned int input = 0; input < 27; input++) {
            cout << dec << setw(2) << DDT[output][input] << " "; 
        }
        cout << endl; 
    }
}

void Sbox::displayOutputCompatibleWith(const Tryte &input) const
{
    vector<TryteSTCompatible>::const_iterator it; 
    
    cout << "Output differences compatible with the input difference "; 
    cout << input << " and the weight of the transition." << endl; 
    for (it = outputDiff[(unsigned int)input].begin(); 
         it != outputDiff[(unsigned int)input].end(); ++it) {
             cout << *it << endl; 
    }
}

void Sbox::displayInputCompatibleWith(const Tryte &output) const
{
    vector<TryteSTCompatible>::const_iterator it; 
    
    cout << "Input differences compatible with the output difference "; 
    cout << output << " and the weight of the transition." << endl; 
    for (it = inputDiff[(unsigned int)output].begin(); 
         it != inputDiff[(unsigned int)output].end(); ++it) {
             cout << *it << endl; 
    }
}

void Sbox::DDTInitialization()
{ 
    // DDT[output][input] = # (x | S(x + input) - S(x) = output) )

    // Initialize the DDT with 0
    for (unsigned int i = 0; i < 27; i++) {
        for (unsigned int j = 0; j < 27; j++)
            DDT[i][j] = 0; 
    }
   
    // Initialize the DDT
    for (unsigned int x = 0; x < 27; x++) { 
        for (unsigned int input = 0; input < 27; input++) {
            Tryte xPlusInput = Tryte(x) + Tryte(input);
            Tryte outputDifference = Tryte(SBOX[(unsigned int)xPlusInput]) - Tryte(SBOX[x]); 
            DDT[(unsigned int)outputDifference][input] += 1;
        }
    }
}

void Sbox::outputInputDiffInitialization()
{
    for (unsigned int output = 0; output < 27; output++) {
        for (unsigned int input = 0; input < 27; input++) {
            if (DDT[output][input] == 3) { 
                outputDiff[input].push_back(TryteSTCompatible(output, Weight(2)));
                inputDiff[output].push_back(TryteSTCompatible(input,  Weight(2)));  
            }
            if (DDT[output][input] == 2) { 
                outputDiff[input].push_back(TryteSTCompatible(output, Weight(0, 1)));
                inputDiff[output].push_back(TryteSTCompatible(input,  Weight(0, 1)));  
            }
            if (DDT[output][input] == 1) {
                outputDiff[input].push_back(TryteSTCompatible(output, Weight(3)));
                inputDiff[output].push_back(TryteSTCompatible(input,  Weight(3)));  
            }
            if (DDT[output][input] == 27) { 
                outputDiff[input].push_back(TryteSTCompatible(output, Weight(0))); 
                inputDiff[output].push_back(TryteSTCompatible(input,  Weight(0))); 
            }
        }
    }
}

bool Sbox::isInKernel(Tryte aTryte1, Tryte aTryte2, Tryte aTryte3)
{
    for (unsigned int i = 0; i < 3; i++) 
        if ((aTryte1.trit(i) + aTryte2.trit(i) + aTryte3.trit(i))% 3 != 0)
            return false; 
    return true;
}

vector<TryteColumnSTCompatible> Sbox::getTryteColumnsBeforeST(unsigned i, 
                                                              unsigned j, 
                                                              unsigned k)
{
    vector<TryteSTCompatible>::iterator itB0; 
    vector<TryteSTCompatible>::iterator itB1;
    vector<TryteSTCompatible>::iterator itB2;
    vector<TryteColumnSTCompatible> res; 
    for (itB0 = inputDiff[i].begin(); itB0 != inputDiff[i].end(); itB0++) {
        for (itB1 = inputDiff[j].begin(); itB1 != inputDiff[j].end(); itB1++) {
            for (itB2 = inputDiff[k].begin(); itB2 != inputDiff[k].end(); itB2++) {
                // we only want the in-kernel box-columns
                if (isInKernel(*itB0, *itB1, *itB2)) { 
                    res.push_back(TryteColumnSTCompatible(*itB0, *itB1, *itB2));
                }
            }
        }
    }
    return res; 
}

void Sbox::inKernelTryteColumnBeforeSTInitialization()
{
    Sbox::inKernelTryteColumnBeforeST.resize(27);
    for (unsigned int i = 0; i < 27; i++) {
        inKernelTryteColumnBeforeST[i] = vector<vector<vector<TryteColumnSTCompatible> > > (i + 1);
        for (unsigned int j = 0; j <= i; j++) {
            inKernelTryteColumnBeforeST[i][j] = vector<vector<TryteColumnSTCompatible> >(j + 1); 
        }
    }

    for (unsigned int i = 0; i < 27; i++) {
        for (unsigned int j = 0; j <= i; j++) {
            for (unsigned int k = 0; k <= j; k++) {
                inKernelTryteColumnBeforeST[i][j][k] = getTryteColumnsBeforeST(i, j, k);
                // The in-kernel box-columns are sorted according to their cost
                sort(inKernelTryteColumnBeforeST[i][j][k].begin(), inKernelTryteColumnBeforeST[i][j][k].end() ); 
            }
        }
    }
}

bool Sbox::areSTCompatible(Tryte inputTryte, Tryte outputTryte, 
                           Weight &accumulatedWeight) const
{
    vector<TryteSTCompatible>::const_iterator it; 
    for (it  = outputDiff[(unsigned int)inputTryte].begin(); 
         it != outputDiff[(unsigned int)inputTryte].end(); it++) {
        if (*it == outputTryte) {
            accumulatedWeight += it->weight; 
            return true; 
        }
    }
    return false; 
}
