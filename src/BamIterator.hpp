// BamIterator.hpp:
// Virtual iterator class with concrete subclasses of
// BamRangeIterator and BamFileIterator.

#ifndef BAMITERATOR_H
#define BAMITERATOR_H

#include "Template.hpp"


class BamIterator {

public:

    const bam_index_t *bindex;
    bam1_t *bam;
    bool iter_done;

    typedef map<string, Template> Templates;
    Templates templates;
    queue<list<const bam1_t *> > complete;
    list<const bam1_t *> invalid;

    // constructor / destructor
    BamIterator(const bam_index_t *bindex) :
        bindex(bindex), bam(NULL), iter_done(false) {}

    virtual ~BamIterator() {
        if (NULL != bam)
            bam_destroy1(bam);
    }

    // process
    void process(const bam1_t *bam) {
        // FIXME: combination of RG and qname?
        const string s = bam1_qname(bam);
        bool mates = templates[s].add_segment(bam, complete, invalid);
        if (mates) {
            if (templates[s].size() == 0)
                templates.erase(s);
        }
    }

    // yield
    void yield(bamFile bfile, bam_mates_t *result) {
        if (complete.empty())
            iterate_inprogress(bfile);
        if (complete.empty())
            finalize_inprogress(bfile);

        list<const bam1_t *> elts;
        bool mated = false;
        if (!complete.empty()) {
            elts = complete.front();
            complete.pop();
            mated = true;
        } else if (!invalid.empty()) {
            elts.push_back(invalid.front());
            invalid.pop_front();
        }

        bam_mates_realloc(result, elts.size(), mated);
        for (int i = 0; !elts.empty(); ++i) {
            result->bams[i] = elts.front();
            elts.pop_front();
        }
    }

    virtual void iterate_inprogress(bamFile bfile) = 0;

    virtual void finalize_inprogress(bamFile bfile) {
        Templates::iterator it;
        // push 'inprogress' to 'invalid'
        for (it = templates.begin(); it != templates.end(); ++it)
            it->second.cleanup(invalid);
        templates.clear();
    }

};

#endif
