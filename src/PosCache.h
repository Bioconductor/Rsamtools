#ifndef POS_CACHE_H
#define POS_CACHE_H

#include <map>
#include <utility>
#include <stdio.h>
#include <iostream>
using std::ostream;
using std::cout;
using std::endl;
using std::map;
using std::pair;
using std::make_pair;

// technically could embed strand information in a nucleotide node,
// but this would hurt readability (assuming only C++98 facilities)
// since the resulting map would require keys of <nucleotide,strand>;
// there's already enough pair.first, pair.second, and make_pair hanky
// panky.
struct StrandTree {
    map<char, int> strandMap;
    StrandTree() : strandMap() { }
    void increment(char strand) {
        // increment count at strand
        strandMap.insert(make_pair(strand, 0)).first->second++;
        //print();
    }
    int numStrands() const { return strandMap.size(); }
    int count() const {
        int count = 0;
        for(map<char,int>::const_iterator it = strandMap.begin();
            it != strandMap.end(); ++it) {
            count += it->second;
        }
        return count;
    }
    // void print() const {
    //     typedef map<char, int>::const_iterator iter;
    //     for(iter it = strandMap.begin(); it != strandMap.end(); ++it) {
    //         printf("strand %c count: %u - ", it->first, it->second);
    //     } printf("\n");
    // }
};

struct PosCache {
    map<char, StrandTree> nucMap;
    PosCache() : nucMap() { }
    // void print() const {
    //     typedef map<char, StrandTree>::const_iterator iter;
    //     for(iter it = nucMap.begin(); it != nucMap.end(); ++it) {
    //         printf("* nuc %c\n", it->first);
    //         it->second.print();
    //     }
    // }
    void increment(char strand, char nucleotide) {
        pair<char, StrandTree> tempNode = make_pair(nucleotide, StrandTree());
        // get tree for nucleotide
        StrandTree &st = nucMap.insert(tempNode).first->second;
        st.increment(strand);
        //print(); printf("\n");
    }
    int numNucs() { return nucMap.size(); }
    int count()  const {
        int count = 0;
        for(map<char,StrandTree>::const_iterator it = nucMap.begin();
            it != nucMap.end(); ++it) {
            count += it->second.count();
        }
        return count;
    }
    int primaryNucCount() const {
        int max = 0;
        for(map<char,StrandTree>::const_iterator it = nucMap.begin();
            it != nucMap.end(); ++it) {
            //printf("nuc %c count %d\n", it->first, it->second.count());
            if(it->second.count() > max)
                max = it->second.count();
        }
        return max;
    }
    void clear() {
        nucMap.clear();
    }
};

#endif // POS_CACHE_H
