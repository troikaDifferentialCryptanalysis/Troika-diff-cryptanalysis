/** trailCore.cpp */ 
#include "trailCore.h"
#include "state.h"
#include "forwardExtension.h"
#include "backwardExtension.h"

TrailCore::TrailCore(const vector<TroikaState>& differences, Weight wMinRev, 
                     Weight wMinDir, const vector<Weight> weights): 
nrRounds(differences.size() / 2 + 1), differences(differences),
wMinRev(wMinRev), wMinDir(wMinDir), weights(weights) 
{
    assert(differences.size() % 2 == 0); 
    assert(weights.size() == nrRounds - 2); 
    weight = wMinRev + wMinDir;
    for (vector<Weight>::const_iterator it = weights.begin(); it != weights.end(); it++)
        weight += (*it); 
}

TrailCore::TrailCore(const TroikaState& stateA, const TroikaState& stateB, 
Weight wMinRevA, Weight wMinDirB): 
nrRounds(2), wMinRev(wMinRevA), wMinDir(wMinDirB)
{
    differences.push_back(stateA);
    differences.push_back(stateB); 
    weight = wMinRevA + wMinDirB;
}

TrailCore::TrailCore(const TroikaState& stateA, const TroikaState& stateB, 
                     const TroikaState& stateC, const TroikaState& stateD, 
                     Weight wMinRevA, Weight wBC, Weight wMinDirD):
    nrRounds(3), wMinRev(wMinRevA), wMinDir(wMinDirD)
{
    differences.push_back(stateA); 
    differences.push_back(stateB); 
    differences.push_back(stateC); 
    differences.push_back(stateD);
    weights.push_back(wBC); 
    weight = wMinRevA + wBC + wMinDirD;
}
        
TrailCore::TrailCore(istream& fin)
{
    load(fin);
}

void TrailCore::clear()
{
    nrRounds = 0; 
    differences.clear(); 
    wMinRev = Weight(0); 
    wMinDir = Weight(0); 
    weights.clear(); 
    weight = Weight(0); 

} 

void TrailCore::load(istream &fin)
{ 
    fin >> hex;  
    fin >> nrRounds;
    if (nrRounds < 2 || fin.eof()) { 
        throw Exception();
    }
    fin >> wMinRev.integer; 
    fin >> wMinRev.logPart;
    weight = wMinRev;

    weights.resize(nrRounds - 2);
    for (unsigned int i = 0; i < nrRounds - 2; i++) { 
        fin >> weights[i].integer;
        fin >> weights[i].logPart; 
        weight += weights[i]; 
    }
    fin >> wMinDir.integer; 
    fin >> wMinDir.logPart; 
    weight += wMinDir; 

    differences.resize(2*(nrRounds - 1));
    for (unsigned int i = 0; i < 2*(nrRounds - 1); i++)
        differences[i].load(fin); 
}

void TrailCore::save(ostream &fout) const
{
    vector<TroikaState>::const_iterator it; 
    fout << hex;
    fout << nrRounds << " "; 
    fout << wMinRev.integer << " ";  
    fout << wMinRev.logPart << " ";  
    for (unsigned int i = 0; i < weights.size(); i++) {
        fout << weights[i].integer << " "; 
        fout << weights[i].logPart << " "; 
    }
    fout << wMinDir.integer << " "; 
    fout << wMinDir.logPart << " "; 
    for (it = differences.begin(); it != differences.end(); it++) 
        (*it).save(fout);
    fout << endl; 
}

void TrailCore::extendForward(const ForwardExtension &extension)
{
    if (nrRounds == 0)
        throw Exception("The trail core is empty. It cannot be extended");

    nrRounds += 1; 
    differences.push_back(extension.stateC); 
    differences.push_back(extension.stateD);
    weights.push_back(extension.wBC);
    weight = weight - wMinDir + extension.wBC + extension.wMinDirD;
    wMinDir = extension.wMinDirD;
}

void TrailCore::extendBackward(const BackwardExtension &extension)
{
    if (nrRounds == 0)
        throw Exception("The trail core is empty. It cannot be extended");

    nrRounds += 1; 
    differences.insert(differences.begin(), extension.stateB); 
    differences.insert(differences.begin(), extension.stateA);
    weights.insert(weights.begin(), extension.wBC);
    weight = extension.wMinRevA + extension.wBC - wMinRev + weight; 
    wMinRev = extension.wMinRevA; 
}

void TrailCore::checkTrailCore(const Sbox &sbox) const
{ 
    // Check the size of the vectors differences and weights 
    if (differences.size() != 2 * (nrRounds - 1)) {
        cerr  << "The size of the vector differences is incorrect." << endl;
        throw Exception("The size of the vector differences is incorrect");

    } 
    if (weights.size() != nrRounds - 2) {
        cerr << "The size of the vector weights is incorrect." << endl;
        throw Exception("The size of the vector weights is incorrect");
    }
    if (nrRounds < 2)
        // the trail is empty
        return; 
    // Check weights 
    if (wMinRev != 2 * differences[0].getNrActiveTrytes()) {
        cerr << "The minimum reverse weight is incorrect; it should be " << 2 * differences[0].getNrActiveTrytes() << " instead of " << wMinRev << "."<< endl; 
        cerr << endl <<  *this  << endl; 
        throw Exception("The minimum reverse weight is incorrect");
    }
    if (wMinDir != 2 * differences.back().getNrActiveTrytes()) {
        cerr << "The minimum weight is incorrect; it should be " << 2 * differences.back().getNrActiveTrytes() << " instead of " << wMinRev << "." << endl; 
        throw Exception("The minimum weight is incorrect");
    }
    // Check compatibility between consecutive states through Lambda 
    TroikaState afterL; 
    for (unsigned int i = 0; i < nrRounds - 1; i++) {
        afterL = differences[2 * i]; 
        afterL.L();
        if (afterL != differences[2 * i + 1]) {
            cerr << "The difference at index " << dec << 2 * i << " is incompatible through Lambda with the difference at index " << dec << 2 * i  + 1 << "." << endl;
            throw Exception("Incompatible states found in the trail core.");
        }
    }
    // Check compatibility between consecutive states through 
    // SubTrytes and the weight of the transition 
    Weight transitionWeight;
    Weight sum = wMinRev + wMinDir; 
    for (unsigned int i = 0; i < nrRounds - 2; i++) {
        if (!sbox.areSTCompatible(differences[ 2 * i + 1 ], differences[2 *(i + 1)], transitionWeight)) {
            cerr << "The difference at index " << dec << 2 * i + 1 << " is incompatible through ST with the difference at index " << dec << 2 * (i  + 1) << "." << endl;
            throw Exception("Incompatible states found in the trail core.");
        }
        if ( transitionWeight != weights[i] ) {
            cerr << "The weight at index " << dec << i  << " is incorrect." << endl;
            throw Exception("Incorrect weight.");

        }
        sum += transitionWeight;
    }
    if (sum != weight) {
        cerr << "The total weight of the trail is incorrect; it should be " << sum << "." << endl;
        throw Exception("The total weight in the trail is incorrect.");

    }
}

bool TrailCore::operator < (const TrailCore& other) const
{
    assert( nrRounds >= 2 && other.nrRounds >= 2 );
    for (unsigned int i = 0; (i < nrRounds - 1) && (i < other.nrRounds - 1); i++) {
        if (differences[2 * i] < other.differences[2 * i])
            return true; 
        else if (other.differences[2 * i] < differences[2 * i])
            return false;  
    }
    return false; 
}

void TrailCore::makeCanonical()
{
    TrailCore currentMin(*this);

    unsigned int dzMin = 0; 

    for (unsigned int dz = 1; dz < SLICES; dz++) {
        TrailCore current(*this); 
        current.translate(dz);
        if (current < currentMin) {
            currentMin = current; 
            dzMin = dz; 
        }
    }
    translate(dzMin);
}

void TrailCore::translate(unsigned int dz)
{
    for (vector<TroikaState>::iterator it = differences.begin(); 
         it != differences.end(); it++) {
        (*it).translate(dz);
    }
}


ostream & operator << (ostream &fout, const TrailCore &aTrailCore)
{
    fout << "A " << aTrailCore.nrRounds ; 
    fout << "-rounds trail core of weight " << aTrailCore.weight << endl << endl;
    fout << "a0 - weight : " << aTrailCore.wMinRev ; 
    fout << aTrailCore.differences[0] << endl; 
    for (unsigned int i = 1; i < aTrailCore.nrRounds - 1; i++) {
        fout << "w (a" << i << " --ST--> b" << i << " ) =  " << aTrailCore.weights[i-1] << endl; 
        fout << "a" << i ; 
        fout << aTrailCore.differences[2*i - 1] << endl;
        fout << "b" << i; 
        fout << aTrailCore.differences[2*i] << endl; 
    }
    fout << "b" << aTrailCore.nrRounds - 1 << " - weight : " << aTrailCore.wMinDir; 
    fout << aTrailCore.differences.back() << endl; 
    return fout;
}

TrailFileIterator::TrailFileIterator(const string& fileName):fin(fileName)
{
    
    if (!fin) 
        throw Exception((string)"File '" + fileName + (string)"' cannot be read." ); 
    end = false; 
    ++(*this);  
}

void TrailFileIterator::operator++()
{
    try {
        current.load(fin);
    } catch (Exception) {
        end = true; 
    }   
}

const TrailCore& TrailFileIterator::operator*() const
{
    return current; 
}

bool TrailFileIterator::isEnd() const
{
    return end; 
}

void produceHumanReadableFile(const string& fileName,
                              bool verbose, 
                              double time)
{
    
    ifstream fin(fileName.c_str());
    string fileName2 = fileName+".txt";
    ofstream fout(fileName2.c_str());
    if (verbose)
        cout << "Writing " << fileName2 << flush;
    
    vector<UINT64> countPerWeight;
    vector<UINT64> countPerLength;
    UINT64 totalCount = 0;
    unsigned int minWeight;
    Sbox sbox;

    if (time)
        fout << "Execution time : " << time << " seconds." << endl; 
   
    while (!(fin.eof())) {
        try {
            TrailCore trail(fin);
            trail.checkTrailCore(sbox);
            if (ceil((double long)trail.weight) >= countPerWeight.size())
                countPerWeight.resize(ceil((double long)trail.weight + 1), 0);
            countPerWeight[ceil((double long)trail.weight)]++;
            if (trail.nrRounds >= countPerLength.size())
                    countPerLength.resize(trail.nrRounds + 1, 0);
                countPerLength[trail.nrRounds]++;
            totalCount++;
        }
        catch(Exception e) {
            cout <<  e.reason << endl; 
        }
    }
    
    if (totalCount == 0) {
        fout << "No trails found in file " << fileName << "!" << endl;
        
    }
    minWeight = 0;
    while ((minWeight < countPerWeight.size()) && (countPerWeight[minWeight] == 0))
        minWeight++;
    for (unsigned int i = 0; i < countPerLength.size(); i++) {
        if (countPerLength[i] > 0)
            fout << dec << countPerLength[i] << " trails of length " << dec << i << " read and checked." << endl;
    }
    fout << "Minimum weight: " << dec << minWeight << endl;
    for (unsigned int i = minWeight; i < countPerWeight.size(); i++) {
        if (countPerWeight[i] > 0) {
            fout.width(8); fout.fill(' ');
            fout << dec << countPerWeight[i] << " trails of weight ]";
            fout.width(2); fout.fill(' ');
            fout << i - 1 << " , " ;
            fout.width(2); fout.fill(' ');
            fout << i << "]" << endl;
        }
    }
    fout << endl; 
}

void produceHumanReadableFileTwoRoundTrailCores(const string& fileName, 
                                                unsigned int alpha, 
                                                unsigned int beta, 
                                                bool verbose, 
                                                double time)
{
    ifstream fin(fileName.c_str());
    string fileName2 = fileName+".txt";
    ofstream fout(fileName2.c_str());
    if (verbose)
        cout << "Writing " << fileName2 << flush;
    
    vector<UINT64> countPerCost;
    UINT64 totalCount = 0;
    unsigned int minCost = 0;
    Sbox sbox;

    if (time)
        fout << "Execution time : " << time << " seconds." << endl; 
    while (!(fin.eof())) {
        try {
            TrailCore trail(fin);
            assert(trail.nrRounds == 2);
            trail.checkTrailCore(sbox);
            unsigned int cost = alpha * (long double)trail.wMinRev 
                               + beta * (long double)trail.wMinDir; 
            if (cost >= countPerCost.size())
                countPerCost.resize( cost + 1, 0);
            countPerCost[cost]++; 
            totalCount++;
        }
        catch(Exception e) {
            cout << e.reason << endl; 
        }
    }
    
    if (totalCount == 0) {
        fout << "No trails found in file " << fileName << "!" << endl;
    }
    minCost = 0;
    while ((minCost < countPerCost.size()) && (countPerCost[minCost] == 0))
        minCost++;

    fout << dec << totalCount << " trails of length 2 read and checked." << endl;
    
    fout << "Minimum cost: " << dec << minCost << endl;
    for (unsigned int i = minCost; i < countPerCost.size(); i++) {
        if (countPerCost[i] > 0) {
            fout.width(8); fout.fill(' ');
            fout << dec << countPerCost[i] << " trails of cost ";
            fout.width(2); fout.fill(' ');
            fout << i << endl ;
        }
    }
    fout << endl;
}
