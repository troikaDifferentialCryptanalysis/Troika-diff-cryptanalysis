/** extensionsIterator.h */ 
#ifndef EXTENSIONS_ITERATOR_H
#define EXTENSIONS_ITERATOR_H

/*

This code uses functions and ideas from KeccakTools
(https://github.com/KeccakTeam/KeccakTools/tree/master/Sources)
and XooTools (https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Sources). 
We thank Silvia Mella and Gilles Van Assche for their intelligible code. 

We define here the general form of the iterator used to extend trail cores 
backward inside the kernel, forward outside the kernel and backward outside 
the kernel (see Appendix C). 
*/

#include "state.h"


template <class Preparation,
          class Parts,  
          class CachedRepresentation,
          class CostFunction,
          class Extension>
          class GenericExtensionIterator
{
    public: 
        /** The state that defines an extension is chosen step by step
          * (see Appendix C). This attribute is a reference to a vector, 
          * initialized by a Preparation, that stores all the parts of the state
          * that must be chosen (eg the trits of the slices, the box-columns). 
          * Each element of the vector exists in different versions. The goal of
          * the iterator is to iterate over the different choices possible for 
          * each element.
          */ 
        vector<Parts>& partsList;
        /** The cache representation of the extension. */ 
        CachedRepresentation cache; 
        /** The cost function.*/ 
        CostFunction costF; 
        /** The output representation of the extension. */ 
        Extension out; 
        /** The maximum weight of the extension. */ 
        long double maxWeightExtension;
        /** Attribute that indicates whether the iterator has reached the end or not. */
        bool end;
        /** Index of the current part of the extension that has to be chosen. */ 
        int indCurPart;
        /** Index of the last part of the extension that has to be chosen,
          * ie partsList.size() - 1. If indCurPart == indLastPart, that means 
          * that the extension is completely chosen.
          */ 
        int indLastPart;
    public:
    /** The constructor.
       * @param  prep   The preparation of the extension that computes all the 
       *                information needed for the extension, including the 
       *                vector partsList.
       * @param  state  If we want to extend a trail core (a1, b1, ... a_{k-1}, b_{k-1})
       *                backward, @state is a1. If we want to extend it forward, @state
       *                is b_{k-1}.
       */
        GenericExtensionIterator(Preparation& prep, const TroikaState& state): 
            partsList(prep.getPartsList()), 
            cache(prep, state), 
            costF(prep, state),
            out(prep, state),
            indCurPart(-1), 
            indLastPart(-1), 
            maxWeightExtension(prep.maxWeightExtension)
        {
            end = ! prep.couldBeExtended();
            if (end == false) {
                indLastPart = partsList.size() - 1;
                indCurPart = -1; 
                ++(*this); 
            }                
        }

	    /** It indicates if there are no more extensions to reach. */
        bool isEnd() const 
        {
            return end; 
        }

        /** It goes (if possible) to the next leaf of the tree that corresponds 
          * to a valid extension (see Appendix C.)
          */
        void operator++()
        {
            if (end)
                return; 
            do {
                do {
                    if (!next()) {
                        end = true; 
                        return; 
                    }
                } while (indCurPart != indLastPart); 
                out.set(partsList, cache, costF);
            } while (out.isValidAndBelowWeight(maxWeightExtension) == false);
        }

        /** It returns a constant reference to the current extension. */
        const Extension& operator*()
        {
            return out; 
        }
    private:
        /** This functions depends on the way we want to prune the tree. 
          * It moves if possible to the next valid node of the tree.
          */
        bool next(); 

        /** It moves if possible to the first child of the current node of 
          * the tree by choosing a value for the next element of partsList
          * (if it exists).
          * @return true if a chilf is found, false otherwise. 
          */
        bool toChild()
        { 
            if (indCurPart == indLastPart)
                return false; 
            indCurPart++; 
            partsList[indCurPart].setFirstValue(cache);
            cache.push(partsList[indCurPart]);
            return true; 
        }
        /** It moves to the parent of the current node by removing the value of 
          * the current element of partsList.
          * @return true if the node has a parent, false if it is the root. 
          */
        bool toParent()
        {
            if (indCurPart == 0)
                return false; 
            cache.pop(partsList[indCurPart]);
            indCurPart --; 
            return true;
        }

        /** It moves to the first sibling of the current node of the tree by 
          * changing the value of the current element of partsList.
          * @return true if a sibling is found, false otherwise.
          */
        bool toSibling()
        { 
            cache.pop(partsList[indCurPart]); 
            if (! partsList[indCurPart].setNextValue(cache) ) {
                cache.push(partsList[indCurPart]); 
                return false; 
            }
            cache.push(partsList[indCurPart]); 
            return true;           
        }
};

#endif
