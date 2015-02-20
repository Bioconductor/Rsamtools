// BamIterator.h:
// Virtual iterator class with concrete subclasses of
// BamRangeIterator and BamFileIterator.

#ifndef BAMITERATOR_H
#define BAMITERATOR_H

#include <set>
#include <map>
#include "Template.h"

class BamIterator {

public:

    const bam_index_t *bindex;
    bam1_t *bam;
    bool iter_done;
    char qname_prefix, qname_suffix;
    bam_qname_f qname_trim;

    typedef map<string, Template> Templates;
    Templates templates;
    queue<list<const bam1_t *> > complete;
    queue<list<const bam1_t *> > ambiguous;
    queue<list<const bam1_t *> > unmated;
    set<string> touched_templates;

    // constructor / destructor
    BamIterator(const bam_index_t *bindex, char qname_prefix, 
        char qname_suffix, bam_qname_f qname_trim) :
        bindex(bindex), bam(NULL), iter_done(false), 
        qname_prefix(qname_prefix), qname_suffix(qname_suffix),
        qname_trim(qname_trim) {}

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
        char *s = qname_trim(bam, qname_prefix, qname_suffix);
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
        // transfer Template::ambiguous to BamIterator::ambiguous
        // transfer Template::inprogress and Template::invalid to 
        // BamIterator::unmated
        for (it = templates.begin(); it != templates.end(); ++it)
            it->second.cleanup(ambiguous, unmated);
        templates.clear();
    }

};

#endif
