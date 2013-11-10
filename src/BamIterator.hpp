// BamIterator.hpp:
// Virtual iterator class with concrete subclasses of
// BamRangeIterator and BamFileIterator.

#ifndef BAMITERATOR_H
#define BAMITERATOR_H

#include <set>
#include "Template.hpp"


class BamIterator {

public:

    const bam_index_t *bindex;
    bam1_t *bam;
    bool iter_done;

    typedef map<string, Template> Templates;
    Templates templates;
    queue<list<const bam1_t *> > complete;
    queue<list<const bam1_t *> > ambiguous;
    queue<list<const bam1_t *> > unmated;
    set<string> touched_templates;

    // constructor / destructor
    BamIterator(const bam_index_t *bindex) :
        bindex(bindex), bam(NULL), iter_done(false) {}

    virtual ~BamIterator() {
        if (NULL != bam)
            bam_destroy1(bam);
    }

    void mate_touched_templates() {
        for (set<string>::iterator it=touched_templates.begin();
             it != touched_templates.end(); ++it) {
            templates[*it].mate(complete);
            if (templates[*it].empty())
                templates.erase(*it);
        }
        touched_templates.clear();
    }

    // process
    void process(const bam1_t *bam) {
        // FIXME: combination of RG and qname?
        const string s = bam1_qname(bam);
        if (templates[s].add_segment(bam))
            touched_templates.insert(s);
    }

    // yield
    void yield(bamFile bfile, bam_mates_t *result) {
        if (complete.empty())
            iterate_inprogress(bfile);
        if (complete.empty())
            finalize_inprogress(bfile);

        list<const bam1_t *> elts;
        int mated = 1;
        if (!complete.empty()) {
            elts = complete.front();
            complete.pop();
        } else if (!ambiguous.empty()) {
            elts = ambiguous.front();
            ambiguous.pop();
            mated = 2;
        } else if (!unmated.empty()) {
            elts = unmated.front();
            unmated.pop();
            mated = 3;
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
        // push 'Template::ambiguous' to 'ambiguous',
        // 'Template::inprogress' and 'Template::invalid' to 'unmated'
        for (it = templates.begin(); it != templates.end(); ++it)
            it->second.cleanup(ambiguous, unmated);
        templates.clear();
    }

};

#endif
