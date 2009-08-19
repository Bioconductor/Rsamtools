#ifndef _RSNAP_H_
#define _RSNAP_H_

#include <Rdefines.h>

typedef struct _snap_t *_SNAP_PTR;

_SNAP_PTR _snap_new();
void _snap_append(_SNAP_PTR ptr, const char *string);
SEXP _snap_as_XStringSet(_SNAP_PTR ptr, const char *baseclass);

#endif
