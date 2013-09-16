// BamIterator.hpp:
// Virtual iterator class with concrete subclasses of
// BamRangeIterator and BamFileIterator.

#ifndef BAMITERATOR_H
#define BAMITERATOR_H

#include <map>
#include <queue>
#include <string>
#include "Template.hpp"
#include "scan_bam_data.h" // new header


class BamIterator {

public:

    const bam_index_t *bindex;
    bam1_t *bam;
    bool iter_done;

    typedef map<string, Template> Templates;
    Templates templates;
    queue<list<const bam1_t *> > complete, incomplete;

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
        bool mates = templates[s].add_segment(bam);
        if (mates) {
            complete.push(templates[s].get_complete());
            if (templates[s].size() == 0)
                templates.erase(s);
        }
    }

    // yield
    void yield(bamFile bfile, bam_mates_t *result) {
        if (complete.size() == 0)
            iterate_complete(bfile);
        if (complete.size() == 0)
            iterate_incomplete(bfile);

        list<const bam1_t *> elts;
        bool mated = false;
        if (complete.size() != 0) {
            elts = complete.front();
            complete.pop();
            mated = true;
        } else if (incomplete.size() != 0) {
            elts = incomplete.front();
            incomplete.pop();
        }

        bam_mates_realloc(result, elts.size(), mated);
        for (int i = 0; !elts.empty(); ++i) {
            result->bams[i] = elts.front();
            elts.pop_front();
        }
    }

    virtual void iterate_complete(bamFile bfile) = 0;

    virtual void iterate_incomplete(bamFile bfile) {
        Templates::iterator it;
        for (it = templates.begin(); it != templates.end(); ++it) {
            // push 'incomplete'
            list<const bam1_t *> tmpl = it->second.cleanup();
            if (tmpl.size())
                incomplete.push(tmpl);
        }
        templates.clear();
    }

};

#endif
