// BamRangeIterator.h:
// Iterator used when reading ranges from a bam file.

#ifndef BAMRANGEITERATOR_H
#define BAMRANGEITERATOR_H

#include "BamIterator.h"

class BamRangeIterator : public BamIterator {

    bam_iter_t iter;

    void iterate_inprogress(bamFile bfile) {
	if ((bam1_t *) NULL == bam) {	// first record 
	    bam = bam_init1();
	    if (bam_iter_read(bfile, iter, bam) < 0) {
		iter_done = true;
		return;
	    }
	}

	do {
	    process(bam);
	    if (bam_iter_read(bfile, iter, bam) < 0)
		iter_done = true;
	} while (!iter_done);
        mate_touched_templates();
    }

    static bool cmp_mpos_template_pair(const pair<int32_t, Template *> m1,
                                       const pair<int32_t, Template *> m2)
    {
        return m1.first < m2.first;
    }

    // mate 'inprogress' segments for each Template in Templates
    void finalize_inprogress(bamFile bfile) {
        int64_t pos = bam_tell(bfile);
        const uint32_t *target_len = header->target_len;
        bam1_t *bam = bam_init1();

        // iterators 
        typedef Templates::iterator template_it;
        typedef list<const bam1_t*>::iterator list_it;
        typedef vector<pair<int32_t, Template *> >::iterator pairs_it;
        typedef map<
            pair<int32_t, int>,
            vector<pair<int32_t, Template *> > >::iterator map_it;

        // groups map:
        // 'key' is the bin (distance) of each mpos from 2^14
        // 'value' is a vector of pairs <mpos, pointer to Template>
        map<pair<int32_t, int>,
            vector<pair<int32_t, Template *> > > groups;
        for (template_it it = templates.begin(); it != templates.end(); ++it) {
            list<const bam1_t *> ip = it->second.get_inprogress();
            for (list_it l_it = ip.begin(); l_it != ip.end(); ++l_it) {
                const bam1_t *curr = *l_it;
                int32_t mpos = curr->core.mpos % target_len[curr->core.mtid];
                pair<int32_t, int> key(curr->core.mtid, mpos / 16384);
                pair<int32_t, Template *> value(curr->core.mpos, &it->second);
                groups[key].push_back(value);
            }
        }

        // Iterate over the groups (bins)
        for (map_it g_it = groups.begin(); g_it != groups.end(); ++g_it) {
            // sort vector of pairs by mpos
            int32_t mtid = g_it->first.first;
            vector<pair<int32_t, Template *> > tmpl_v = g_it->second;
            sort(tmpl_v.begin(), tmpl_v.end(), cmp_mpos_template_pair);

            int32_t
                min_mpos = tmpl_v.front().first,
                max_mpos = tmpl_v.back().first + 1;
            bam_iter_t iter = bam_iter_query(bindex, mtid, min_mpos, max_mpos);
            pairs_it p_it = tmpl_v.begin();
            do {
                /* find the next read / template pair */
                int status = 0;
                while ((status = bam_iter_read(bfile, iter, bam)) >= 0)
                    if (bam->core.pos >= p_it->first)
                        break;
                if (status < 0)
                    break;

                while (p_it != tmpl_v.end() && bam->core.pos > p_it->first)
                    ++p_it;
                if (p_it == tmpl_v.end())
                    break;
                if (p_it->first != bam->core.pos)
                    continue;

                /* move p_it to back of valid range */
                /* FIXME: really only needs to be done if p_it changed */
                pairs_it p_start, p_end;
                p_start = p_end = p_it;
                while (p_end != tmpl_v.end() && p_end->first == p_start->first)
                    ++p_end;

                /* find mates for this read, across all mates */
                for (pairs_it p = p_start; p != p_end; ++p) {
                    // Each Template has an 'inprogress' list that may have > 1 
                    // record. In this case the 'groups' map will have
                    // multiple pairs with different mpos but a pointer to
                    // the same Template. All records in 'inprogress' are
                    // processed when the Template is visited the first time.
                    // To avoid revisiting a Template we could try setting
                    // a global such as 'touched = true'.
                    Template *tmpl = p->second;
                    list<const bam1_t *> ip = tmpl->get_inprogress();
                    for (list_it l_it = ip.begin(); l_it != ip.end(); ++l_it) {
                        const bam1_t *curr = *l_it;
                        char pre = qname_prefix_end();
                        char suf = qname_suffix_start();
                        const char *tmpl_q;
                        tmpl_q = tmpl->qname_trim(bam1_qname(curr), pre, suf);
                        const char *mate_q;
                        mate_q = tmpl->qname_trim(bam1_qname(bam), pre, suf);
                        if (tmpl->is_valid(bam) && 
                            tmpl->is_template(tmpl_q, mate_q, bam) &&
                            tmpl->is_mate(curr, bam, target_len)) {
                            tmpl->add_segment(bam);
                            break;
                        }
                    }
                }
            } while (p_it != tmpl_v.end());
            bam_iter_destroy(iter);
        }
        bam_destroy1(bam);

        // mate Templates
        for (template_it it = templates.begin(); it != templates.end(); ++it)
            it->second.mate(complete, target_len);
            
        BamIterator::finalize_inprogress(bfile);
        (void) bam_seek(bfile, pos, SEEK_SET);
    }

public:

    // constructor / destructor
    BamRangeIterator(bamFile bfile, const bam_index_t *bindex,
                     int32_t tid, int32_t beg, int32_t end) :
        BamIterator(bfile, bindex)
    {
	iter = bam_iter_query(bindex, tid, beg, end);
    }

    ~BamRangeIterator() {
	bam_iter_destroy(iter);
   }
};

#endif
