#ifndef HTS_UTILITIES_H
#define HTS_UTILITIES_H

#include <stdint.h>
#include <sys/types.h>
#include <htslib/hts.h>

#ifdef __cplusplus
extern "C" {
#endif

htsFormat _hts_utilities_format(const char *filename);

int _hts_utilities_seek(htsFile *fd, off_t offset, int whence);

int64_t _hts_utilities_tell(htsFile *fd);

#ifdef __cplusplus
}
#endif

#endif
