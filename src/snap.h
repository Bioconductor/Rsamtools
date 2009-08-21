#ifndef _RSNAP_H_
#define _RSNAP_H_

#include <Rdefines.h>

typedef struct _snap_t _SNAP_T;

_SNAP_T *_snap_new();
void _snap_append(_SNAP_T *ptr, const char *string);
SEXP _snap_as_XStringSet(_SNAP_T *ptr, const char *baseclass);
SEXP _snap_as_PhredQuality(_SNAP_T *ptr);
void _snap_delete(_SNAP_T *);

void _snap_cleanup_all(); /* only for error handling */

#endif
