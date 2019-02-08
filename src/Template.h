#ifndef TEMPLATE_H
#define TEMPLATE_H

#include <map>
#include <queue>
#include <string>
#include <algorithm>
#include <list>
#include <sam.h>
#include "scan_bam_data.h"

using namespace std;

class Template {

    typedef list<const bam1_t *> Segments;
    typedef Segments::iterator iterator;
    typedef Segments::const_iterator const_iterator;

    Segments inprogress, ambiguous, invalid;

    void add_to_complete(const bam1_t *bam, const bam1_t *mate,
                         queue<Segments> &complete) {
        Segments tmp;
        if (bam->core.flag & BAM_FREAD1) {
            tmp.push_back(bam);
            tmp.push_back(mate);
        } else {
            tmp.push_back(mate);
            tmp.push_back(bam);
        } 
        complete.push(tmp);
    }

public:

    bool touched;

    // check readgroup and trimmed_qname
    bool is_template(const string trimmed_qname,
                     const string m_trimmed_qname,
                     const bam1_t *mate) const {
        bool test = false;

        /* read group */
        const uint8_t
            *aux = bam_aux_get(inprogress.front(), "RG"),
            *m_aux = bam_aux_get(mate, "RG");
        const char
            *rg = (aux == (const uint8_t *) NULL) ?
                (const char *) NULL : bam_aux2Z(aux),
            *m_rg = (m_aux == (const uint8_t *)  NULL) ?
                (const char *) NULL : bam_aux2Z(m_aux);

        if ((aux == NULL && m_aux == NULL) ||
            (aux != NULL && m_aux != NULL && strcmp(rg, m_rg) == 0))
            test = true;

        /* trimmed qname */
        return test && trimmed_qname.compare(m_trimmed_qname) == 0;
    }
    // is_valid checks the following bit flags:
    // 1. Bit 0x1 (multiple segments) is 1 
    // 2. Bit 0x4 (segment unmapped) is 0
    // 3. Bit 0x8 (next segment unmapped) is 0
    // 4. 'mpos' != -1 (i.e. PNEXT = 0)
    bool is_valid(const bam1_t *bam) {
        const bool multi_seg = bam->core.flag & BAM_FPAIRED;
        const bool seg_unmapped = bam->core.flag & BAM_FUNMAP;
        const bool mate_unmapped = bam->core.flag & BAM_FMUNMAP;

        return  multi_seg && !seg_unmapped && 
                !mate_unmapped && (bam->core.mpos != -1);
    }

    // is_mate checks the following bit flags:
    // 1. Bit 0x40 and 0x80: Segments are a pair of first/last OR
    //    neither segment is marked first/last
    // 2. Bit 0x100: Both segments are secondary OR both not secondary
    // 3. Bit 0x10 and 0x20: Strand flag 0x20 of one mate must match strand
    //                       flag 0x10 of the other mate and vice versa
    // 4. Bit 0x2: Both proper OR both not proper 
    // 5. mpos match:
    //      bit 0x10 of rec1 == bit 0x20 of rec2 AND
    //      bit 0x10 or rec2 == bit 0x20 of rec1
    //      segment2 mpos matches segment1 pos
    // 6. tid match
    bool is_mate(const bam1_t *bam, const bam1_t *mate,
                 const uint32_t *target_len) const {
        const bool bam_read1 = bam->core.flag & BAM_FREAD1;
        const bool mate_read1 = mate->core.flag & BAM_FREAD1;
        const bool bam_read2 = bam->core.flag & BAM_FREAD2;
        const bool mate_read2 = mate->core.flag & BAM_FREAD2;
        const bool bam_secondary = bam->core.flag & BAM_FSECONDARY;
        const bool mate_secondary = mate->core.flag & BAM_FSECONDARY;
        const bool bam_proper = bam->core.flag & BAM_FPROPER_PAIR;
        const bool mate_proper = mate->core.flag & BAM_FPROPER_PAIR;
        const bool bam_rev = bam->core.flag & BAM_FREVERSE;
        const bool mate_rev = mate->core.flag & BAM_FREVERSE;
        const bool bam_mrev = bam->core.flag & BAM_FMREVERSE;
        const bool mate_mrev = mate->core.flag & BAM_FMREVERSE;
        const uint32_t
            pos = bam->core.pos % target_len[bam->core.tid],
            mpos = bam->core.mpos % target_len[bam->core.mtid],
            mate_pos = mate->core.pos % target_len[mate->core.tid],
            mate_mpos = mate->core.mpos % target_len[mate->core.mtid];
        return
            ((bam_read1 ^ bam_read2) && (mate_read1 ^ mate_read2)) &&
            (bam_read1 != mate_read1) &&
            (bam_secondary == mate_secondary) &&
            (((bam_rev != mate_mrev) && (bam_mrev != mate_rev)) || 
            ((bam_rev == mate_mrev) && (bam_mrev == mate_rev))) &&
            (bam_proper == mate_proper) &&
            (pos == mate_mpos) && (mpos == mate_pos) &&
            (bam->core.mtid == mate->core.tid);
    }

    Template() { touched = false; }

    bool empty() const {
        return inprogress.empty() && invalid.empty() && ambiguous.empty();
    }
 
    void touched_true() {
        touched = true;
    }

    bool get_touched() {
        return touched;
    }

    list<const bam1_t*> get_inprogress() {
        return inprogress;
    }

    static const char *qname_trim(char *qname, const char prefix,
                                  const char suffix)
    {
        char *end = qname + strlen(qname);

        if (suffix != '\0')
            for (char *e = end; e >= qname; e--) {
                if (*e == suffix) {
                    *e = '\0';
                    end = e;
                    break;
                }
            }

        if (prefix != '\0')
            for (char *s = qname; *s != '\0'; s++) {
                if (*s == prefix) {
                    memmove(qname, s + 1, end - s);
                    break;
                }
            }

        return qname;
    }

    // Returns true if potential mate, false if invalid
    bool add_segment(const bam1_t *bam1) {
        bam1_t *bam = bam_dup1(bam1);
        if (!is_valid(bam)) {
            invalid.push_back(bam);
            return false;
        }
        inprogress.push_back(bam);
        return true;
    }

    void mate(queue<Segments> &complete, const uint32_t *target_len) {
        const int unmated=-1, multiple=-2, processed=-3;
        vector<pair<int, const bam1_t *> >
            status(inprogress.size(),
                   pair<int, const bam1_t *>(unmated, (const bam1_t *) NULL));
        Segments::iterator it0;

        // identify unambiguous and ambiguous mates
        it0 = inprogress.begin();
        for (unsigned int i = 0; i < inprogress.size(); ++i) {
            status[i].second = *it0;
            Segments::iterator it1 = it0;
            for (unsigned int j = i + 1; j < inprogress.size(); ++j) {
                ++it1;
                if (is_mate(*it0, *it1, target_len)) {
                    status[i].first = status[i].first == unmated ? j : multiple;
                    status[j].first = status[j].first == unmated ? i : multiple;
                }
            }
            ++it0;
        }

        // process unambiguous and ambiguous mates
        for (unsigned int i = 0; i < status.size(); ++i) {
            if (status[i].first == unmated)
                continue;
            if (status[i].first >= 0 && status[status[i].first].first >= 0) {
                // unambiguous mates
                add_to_complete(status[i].second,
                                status[status[i].first].second,
                                complete);
                status[status[i].first].first = processed;
                status[i].first = processed;
            } else if (status[i].first != processed) {
                // ambiguous mates, added to 'ambiguous' queue
                ambiguous.push_back(status[i].second);
                status[i].first = processed;
            }
            ++it0;
        }

        // remove segments that have been assigned to complete or
        // ambiguous queue
        it0 = inprogress.begin();
        for (unsigned int i = 0; i != status.size(); ++i) {
            if (status[i].first == processed) {
                it0 = inprogress.erase(it0);
            } else {
                ++it0;
            }
        }
    }

    // move 'ambiguous' to ambiguous_queue 
    // move 'inprogress' and 'invalid' to 'invalid_queue'
    void cleanup(queue<Segments> &ambiguous_queue,
                 queue<Segments> &invalid_queue) {
        if (!ambiguous.empty())
            ambiguous_queue.push(ambiguous);
        if  (!invalid.empty())
            inprogress.splice(inprogress.end(), invalid);
        if (!inprogress.empty()) {
            invalid_queue.push(inprogress);
            inprogress.clear();
        }
    }

};

#endif
