// BamIterator.h:
// Virtual iterator class with concrete subclasses of
// BamRangeIterator and BamFileIterator.

#ifndef BAMITERATOR_H
#define BAMITERATOR_H

#include <set>
#include <map>
#include "Template.h"
#include "bam_data.h"

class BamIterator {

public:

    const bam_index_t *bindex;
    bam1_t *bam;
    bool iter_done;
    BAM_DATA bam_data;
    char qname_prefix, qname_suffix;

    typedef map<string, Template> Templates;
    Templates templates;
    queue<list<const bam1_t *> > complete;
    queue<list<const bam1_t *> > ambiguous;
    queue<list<const bam1_t *> > unmated;
    set<string> touched_templates;

    // constructor / destructor
    BamIterator(const bam_index_t *bindex, BAM_DATA bam_data) :
        bindex(bindex), bam(NULL), iter_done(false),
        bam_data(bam_data), qname_prefix(bam_data->qnamePrefixEnd),
        qname_suffix(bam_data->qnameSuffixStart) {}

    virtual ~BamIterator() {
        if (NULL != bam)
            bam_destroy1(bam);
    }

    void mate_touched_templates() {
        Rprintf("mate_touched_templates start size: %d\n", templates.size());
        for (set<string>::iterator it=touched_templates.begin();
             it != touched_templates.end(); ++it) {
            templates[*it].mate(complete);
            if (templates[*it].empty())
                templates.erase(*it);
        }
        Rprintf("mate_touched_templates end size: %d\n", templates.size());
        touched_templates.clear();
    }

    // process
    void process(const bam1_t *bam) {
        if (!_filter1_BAM_DATA(bam, bam_data))
            return;
        const char *trimmed_qname =
            Template::qname_trim(bam1_qname(bam), qname_prefix, qname_suffix);
        if (templates[trimmed_qname].add_segment(bam, trimmed_qname))
            touched_templates.insert(trimmed_qname);
    }

    // yield
    void yield(bamFile bfile, bam_mates_t *result) {
        if (complete.empty()) {
            Rprintf("iterate_inprogress\n");
            iterate_inprogress(bfile);
            Rprintf("iterate_inprogress complete.size(): %d\n", complete.size());
        }
        if (complete.empty()) {
            Rprintf("finalize_inprogress\n");
            finalize_inprogress(bfile);
            Rprintf("finalize_inprogress complete.size(): %d\n",
                    complete.size());
        }

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
