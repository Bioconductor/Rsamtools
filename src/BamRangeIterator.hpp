#ifndef BAMRANGEITERATOR_H
#define BAMRANGEITERATOR_H

#include <map>
#include <queue>
#include <string>
#include "Template.hpp"


class BamRangeIterator {

    const bam_index_t *bindex;

    bam_iter_t iter;
    bam1_t *bam;
    bool mated, iter_done;

    typedef map<string, Template> Templates;
    Templates templates;
    queue<Template> pending, incomplete;

    // FPAIRED is also checked in Template (is_mate)
    void process(bam1_t *bam) {
	if (bam->core.flag & BAM_FPAIRED) {
	    const string s = bam1_qname(bam);
	    bool complete = templates[s].add_segment(bam);
	    if (complete) {
		pending.push(templates[s]);
		templates.erase(s);
	    }
	} else {
	    pending.push(Template(bam));
	}
    }
   
    void mate_all(bamFile bfile) {
	int64_t pos = bam_tell(bfile);
	Templates::iterator it = templates.begin();
	while (it != templates.end()) {
	    if (it->second.mate(bfile, bindex)) {
		pending.push(it->second);
		templates.erase(it);
		break;
	    } else {
		incomplete.push(it->second);
		Templates::iterator it0 = it++;
		templates.erase(it0);
	    }
	}
	bam_seek(bfile, pos, SEEK_SET);
    }

    void iterate(bamFile bfile) {
	if (iter_done)
	    return;
	if (NULL == bam) {	// first record 
	    bam = bam_init1();
	    if (bam_iter_read(bfile, iter, bam) < 0) {
		iter_done = true;
		return;
	    }
	}

	int tid = -1, pos = -1;
	bool done = false;
	do {
	    process(bam);
	    bam = bam_init1();
	    if (bam_iter_read(bfile, iter, bam) < 0) {
		iter_done = done = true;
	    } else if (pending.size() != 0) {
		// stop if something to yield AND finished position
		if (tid == -1) { // tid & pos of first pending record
		    tid = bam->core.tid;
		    pos = bam->core.pos;
		} else {
		    done = (bam->core.tid != tid) || (bam->core.pos != pos);
		}
	    }
	} while (!done);
    }

public:

    BamRangeIterator(const bam_index_t *bindex,
		     int tid, int beg, int end, bool mated) :
	bindex(bindex), mated(mated), bam(NULL), iter_done(false)
    {
	iter = bam_iter_query(bindex, tid, beg, end);
    }

    ~BamRangeIterator() {
	bam_iter_destroy(iter);
	if (NULL != bam)
	    bam_destroy1(bam);
    }

    list<bam1_t *> yield(bamFile bfile) {
	if (pending.size() == 0)
	    iterate(bfile);
	if (pending.size() == 0)
	    mate_all(bfile);

	list<bam1_t *> result;
	if (pending.size() != 0) {
	    result = pending.front().values();
	    pending.pop();
	} else if (incomplete.size() != 0) {
	    result = incomplete.front().values();
	    incomplete.pop();
	}

	return result;
    }
};

#endif
