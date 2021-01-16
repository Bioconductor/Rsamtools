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

    BAM_DATA bam_data;

    queue<list<const bam1_t *> > ambiguous;
    queue<list<const bam1_t *> > unmated;
    set<string> touched_templates;

protected:

    typedef map<string, Template> Templates;
    Templates templates;
    queue<list<const bam1_t *> > complete;

    const bam_index_t *bindex;
    bam_header_t *header;
    bam1_t *bam;

    char qname_prefix_end() const {
        if (bam_data == (BAM_DATA) NULL)
            Rf_error("[qname_prefix_end] report to maintainer('Rsamtools')");
        return bam_data->qnamePrefixEnd;
    }

    char qname_suffix_start() const {
        if (bam_data == (BAM_DATA) NULL)
            Rf_error("[qname_suffix_start] report to maintainer('Rsamtools')");
        return bam_data->qnameSuffixStart;
    }

    void mate_touched_templates() {
        for (set<string>::iterator it=touched_templates.begin();
             it != touched_templates.end(); ++it) {
            templates[*it].mate(complete, header->target_len);
            if (templates[*it].empty())
                templates.erase(*it);
        }
        touched_templates.clear();
    }

    // process
    void process(const bam1_t *bam) {
        if (bam_data == (BAM_DATA) NULL)
            Rf_error("[process] report to maintainer('Rsamtools')");
        if (!_filter1_BAM_DATA(bam, bam_data))
            return;
        const char *trimmed_qname =
            Template::qname_trim(bam1_qname(bam), qname_prefix_end(),
                                 qname_suffix_start());
        if (templates[trimmed_qname].add_segment(bam))
            touched_templates.insert(trimmed_qname);
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

public:

    bool iter_done;

    // constructor / destructor
    BamIterator(bamFile bfile, const bam_index_t *bindex) :
        bam_data((BAM_DATA) NULL), bindex(bindex),
        bam((bam1_t *) NULL), iter_done(false)
    {
        (void) bam_seek(bfile, 0, 0);
        header = bam_header_read(bfile);
    }

    virtual ~BamIterator() {
        if ((bam1_t *) NULL != bam)
            bam_destroy1(bam);
        bam_header_destroy(header);
    }

    void set_bam_data(BAM_DATA bd) {
        this->bam_data = bd;
    }

    // yield
    void yield(bamFile bfile, bam_mates_t *result) {
        if (complete.empty() && !iter_done)
            iterate_inprogress(bfile);
        if (complete.empty() && !templates.empty())
            finalize_inprogress(bfile);

        list<const bam1_t *> elts;
        MATE_STATUS mated = MATE_UNKNOWN;
        if (!complete.empty()) {
            elts = complete.front();
            complete.pop();
            mated = MATE_MATED;
        } else if (!ambiguous.empty()) {
            elts = ambiguous.front();
            ambiguous.pop();
            mated = MATE_AMBIGUOUS;
        } else if (!unmated.empty()) {
            elts = unmated.front();
            unmated.pop();
            mated = MATE_UNMATED;
        }

        bam_mates_realloc(result, elts.size(), mated);
        for (int i = 0; !elts.empty(); ++i) {
            result->bams[i] = elts.front();
            elts.pop_front();
        }
    }

};

#endif
