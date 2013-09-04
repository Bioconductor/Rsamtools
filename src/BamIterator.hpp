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
    queue<list<bam1_t *> > complete, incomplete;

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
    void yield(bamFile bfile, bam_mates_t *result, bool do_mate_all) {
        if (complete.size() == 0)
            iterate(bfile);
        if (complete.size() == 0) {
            if (do_mate_all)
                mate_all(bfile);
            else
                append_incomplete();
        } 

        list<bam1_t *> elts;
        if (complete.size() != 0) {
            elts = complete.front();
            complete.pop();
            bam_mates_realloc(result, elts.size());
            result->mates = true;
        } else if (incomplete.size() != 0) {
            elts = incomplete.front();
            incomplete.pop();
            bam_mates_realloc(result, elts.size());
            result->mates = false;
        } else {
            bam_mates_realloc(result, 0);
        }

        for (int i = 0; !elts.empty(); ++i) {
            result->bams[i] = elts.front();
            elts.pop_front();
        }
    }

    // mate_all (BamRangeIterator only)
    void mate_all(bamFile bfile) {
        int64_t pos = bam_tell(bfile);
        Templates::iterator it;
        for (it = templates.begin(); it != templates.end(); ++it) {
            // mate all segments in 'inprogress'
            while (it->second.mate(bfile, bindex)) {
                complete.push(it->second.get_complete());
            }
            // push 'incomplete'
            list<bam1_t *> tmpl = it->second.cleanup();
            if (tmpl.size())
                incomplete.push(tmpl);
        }
        templates.clear();
        bam_seek(bfile, pos, SEEK_SET);
    }

    // append_incomplete (BamFileIterator only)
    void append_incomplete() {
        Templates::iterator it;
        for (it = templates.begin(); it != templates.end(); ++it) {
            list<bam1_t *> tmpl = it->second.cleanup();
            if (tmpl.size())
                incomplete.push(tmpl);
        }
        templates.clear();
    }

    virtual void iterate(bamFile bfile) = 0;
};

#endif
