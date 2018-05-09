#ifndef _COMPAT_BCF_HDR_READ_H_
#define _COMPAT_BCF_HDR_READ_H_

#include <htslib/hts.h>
#include <htslib/vcf.h>

bcf_hdr_t *COMPAT_bcf_hdr_read(htsFile *hfp);

#endif                      /* _COMPAT_BCF_HDR_READ_H_ */
