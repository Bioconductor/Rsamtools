#ifndef TEMPLATE_H
#define TEMPALTE_H

#include <algorithm>
#include <list>
#include "samtools/sam.h"

using namespace std;

class Template {
    typedef list<bam1_t *> Segments;
    typedef Segments::iterator iterator;
    typedef Segments::const_iterator const_iterator;

    char *rg, *qname;
    Segments segments;
    bool complete;

    // FIXED: handle case of rg == rg0 == null
    const int readgroup_q(const bam1_t *mate) const {
        if (rg == '\0')
            return -1;
	uint8_t *aux = bam_aux_get(mate, "RG");
	char *rg0 = '\0';
	if (aux != 0)
	    rg0 = bam_aux2Z(aux);
	return strcmp(rg, rg0);
    }

    const int qname_q(const bam1_t *mate) const {
	return strcmp(qname, bam1_qname(mate));
    }

    bool is_template(const bam1_t *mate) const {
	return (readgroup_q(mate) == 0) && (qname_q(mate) == 0);
    }

    // compare to Herve's is_a_pair() in utilities.c 
    bool is_mate(const bam1_t *bam, const bam1_t *mate) const {
	const bool mate_unmapped = bam->core.flag & BAM_FMUNMAP;
	const bool unmapped = mate->core.flag & BAM_FUNMAP;
	return
	    (bam->core.flag & mate->core.flag & BAM_FPAIRED) &&
	    (bam != mate) &&
	    (mate_unmapped == unmapped) &&
	    (bam->core.mtid == mate->core.tid) &&
	    (bam->core.mpos == mate->core.pos);
    }

    bool needs_mate (const bam1_t *curr) const {
	if (!(curr->core.flag & BAM_FPAIRED))
	    // SINGLE_END
	    return false;

	if (curr->core.mpos == -1)
	    // MATE_UNAVAILABLE
	    return true;
    
	/* scan for previously retrieved */
	for (const_iterator it = segments.begin(); it != segments.end(); ++it) {
	    if (is_mate(curr, *it))
		// MATE_FOUND
		return false;
	}

	return true;
    }

    void find_mate(const bam1_t *curr, bamFile bfile, const bam_index_t * bindex)
    {
	const int tid = curr->core.mtid;
	const int beg = curr->core.mpos;

	if (beg == -1)
	    return;		// no mate available

	bam1_t *bam = bam_init1();
	bam_iter_t iter = bam_iter_query(bindex, tid, beg, beg + 1);
	while (bam_iter_read(bfile, iter, bam) >= 0) {
	    if ((is_template(bam)) && is_mate(curr, bam))
		add_segment(bam_dup1(bam));
}
	bam_destroy1(bam);
	bam_iter_destroy(iter);
    }

public:

    Template() : complete(false), rg('\0') {}

    Template(bam1_t *bam) : complete(true), rg('\0') { add_segment(bam); }

    size_t size() const { return segments.size(); }

    bool add_segment(bam1_t *segment) {
	segments.push_back(segment);
	if (size() == 1) {
	    qname = bam1_qname(segments.front());
	    uint8_t *aux = bam_aux_get(segments.front(), "RG");
	    if (aux != 0)
		rg = bam_aux2Z(aux);
	}
	return is_complete();
    }

    const list<bam1_t *>& values() const { return segments; }

    bool is_complete() {
	if (!complete) {
	    for (iterator it = segments.begin(); it != segments.end(); ++it)
		if (needs_mate(*it))
		    return false;
	    complete = true;
	}
	return complete;
    }

    bool mate(bamFile bfile, const bam_index_t * bindex) {
	if (complete)
	    return complete;
	for (iterator it = segments.begin(); it != segments.end(); ++it)
	    if (needs_mate(*it))
		find_mate(*it, bfile, bindex);
	return is_complete();
    }
};

#endif
