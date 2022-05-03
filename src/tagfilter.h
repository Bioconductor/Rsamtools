#ifndef TAG_FILTER_H
#define TAG_FILTER_H

#include <samtools-1.7-compat.h>
#include <Rdefines.h>

typedef enum { TAGFILT_T_UNSET = 0, TAGFILT_T_INT,
               TAGFILT_T_STRING } TagFilterType;

typedef struct {
    int len;
    TagFilterType type;
    void* ptr;
    
} _TAGFILTER_ELT, *TAGFILTER_ELT;

typedef struct {
    int len;
    const char **tagnames;
    _TAGFILTER_ELT *elts;
} _C_TAGFILTER, *C_TAGFILTER;

C_TAGFILTER _tagFilter_as_C_types(SEXP tl);
void _Free_C_TAGFILTER(C_TAGFILTER ctf);

int _tagfilter(const bam1_t * bam, C_TAGFILTER tagfilter, int irec);

#endif /* TAG_FILTER_H */
