/** trailCore.h*/ 
#ifndef TRAILCORE_H
#define TRAILCORE_H

/*

This code uses functions and ideas from KeccakTools
(https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

The classes of this file are used to store trail cores, to verify their 
consistency, to save them in a file and read them from a file. 
*/

#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include "types.h"
#include "state.h"
#include "sbox.h"

class ForwardExtension;
class BackwardExtension;

/** Class used to store a trail core.
  * A k-round trail core is of the form : 
  * a1 --Lambda--> b1 --ST--> , ..., a_{k-1}--Lambda--> b_{k-1}
  */
class TrailCore
{
    public: 
        /** Number of rounds of the trail core.
          * eg. A--Lambda--> B is a 2-round trail core.
          *     A --Lambda--> B --ST--> C --Lambda--> D is a 3-round trail core. 
          */
        unsigned int nrRounds;
        /** Vector of size 2*(nrRounds - 1) which stores the trail core's 
          * differences. The trail core is 
          * differences[0] --L --> differences[1]--ST--> differences[2]
          * ...  --L-> differences.back()
          */
        vector<TroikaState> differences; 
        /** Minimum reverse weight of the state differences[0] */
        Weight wMinRev;
        /** Minimum weight of the state differences.back() */
        Weight wMinDir;
        /** Vector of size nrRounds - 2 which stores the weights of 
          * the trail core's differentials (all the weights except the
          * minimum reverse weight of differences[0] and the minimum
          * weight of differences.back() ). 
          */
        vector<Weight> weights;
        /** Total weight of the trail core. */
        Weight weight;
    public: 
        TrailCore():nrRounds(0), wMinRev(0), wMinDir(0), weight(0){}
        TrailCore(const vector<TroikaState>& differences, 
                  Weight wMinRev, 
                  Weight wMinDir,
                  const vector<Weight> weights);
        
        TrailCore(const TrailCore& other): nrRounds(other.nrRounds),
         differences(other.differences), wMinRev(other.wMinRev),
         wMinDir(other.wMinDir), weights(other.weights), weight(other.weight){}
        /** Constructor for a 2-round trail core A --Lambda--> B */ 
        TrailCore(const TroikaState& stateA, const TroikaState& stateB, 
                  Weight wMinRevA, Weight wMinDirB);
        /** Constructor for a 3-round trail core
          * A--Lambda--> B --ST--> C --Lambda--> D 
          */ 
        TrailCore(const TroikaState& stateA, const TroikaState& stateB, 
                  const TroikaState& stateC, const TroikaState& stateD, 
                  Weight wMinRevA, Weight wBC, Weight wMinDirD);
        TrailCore(istream& fin);
        /** This constructor initializes a trail core by copying another 
          * trail core given in parameter and extending it with an extension 
          * also given in parameter. 
          * @param trailToExtend The trail core to extend.
          * @param extension     The backward extension used to extend @trailToExtend.
          */
        TrailCore(const TrailCore& trailToExtend, const BackwardExtension& extension)
        :TrailCore(trailToExtend) {extendBackward(extension);}
        /** Same as above but with a forward extension. */ 
        TrailCore(const TrailCore& trailToExtend, const ForwardExtension& extension)
        :TrailCore(trailToExtend) {extendForward(extension);}
        void clear(); 
        void load(istream& fin);
        void save(ostream& fout) const;
        void extendForward(const ForwardExtension& extension);
        void extendBackward(const BackwardExtension& extension);
        /** This method checks the consistency of the trail core.
          * If an inconsistency is found, an Exception is thrown and details
          * about the inconsistency are displayed in the error console (cerr).
          * The aspects tested are:
          * - if the propagation weights declared in the trail match the
          *   propagation weights of the specified differences;
          * - if between two rounds, the specified differences are compatible.
          * @param sbox Sbox used to check if between two rounds, the specified
          *             differences are compatible and the weight of the transition.
          */
        void checkTrailCore(const Sbox& sbox) const;
        /** An arbitrary order relation on trail cores. */ 
        bool operator < (const TrailCore& other) const;
         /* It makes the trail core canonical w.r.t translation along the z-axis.
          */
        void makeCanonical();
        /** It translates the trail core along the z-axis by @dz */ 
        void translate(unsigned int dz);
        friend ostream & operator << (ostream &fout, const TrailCore& aTrailCore);
};

/** Class used to read trail cores from a file. */ 
class TrailFileIterator 
{
    protected: 
        ifstream fin;
        bool end;
        TrailCore current;
    public: 
        /** The constructor of the iterator.
          * @param  aFileName   The name of the file to read from.
          */
        TrailFileIterator(const string& aFileName); 
        /**  It moves the iterator to the next trail in the set. */ 
        void operator++(); 
        /** It gives a constant reference to the current trail. */ 
        const TrailCore& operator*() const;
        /** It indicates whether the iterator has reached the end of the set of trails. */
        bool isEnd() const;
};

/** It reads all the trail cores  in a file, checks their consistency and
  * then produces a report.
  * The report is output in a file with the same file name plus ".txt".
  * @param   fileName   The name of the file containing the trail cores.
  * @param   verbose    If true, the function will display the name of
  *                     the file written to cout.
  * @param   time       A time to write in the report.
  */
void produceHumanReadableFile(const string& fileName,
                              bool verbose = true,
                              double time = 0);

/** It reads all the 2-round trail cores  in a file, checks their consistency 
  * and then produces a report about their costs, of the form 
  * alpha * wMinRev + beta * wMinDir. (see Section 4.2 Generating 2-round trail
  * cores as a tree traversal)
  * The report is output in a file with the same file name plus ".txt".
  * @param   fileName   The name of the file containing the trail cores.
  * @param   alpha      Weighting for the cost. 
  * @param   beta       Weighting for the cost. 
  * @param   verbose    If true, the function will display the name of
  *                     the file written to cout.
  * @param   time       A time to write in the report.
  */
void produceHumanReadableFileTwoRoundTrailCores(const string& fileName, 
                                                unsigned int alpha, 
                                                unsigned int beta, 
                                                bool verbose = true, 
                                                double time = 0);

#endif
