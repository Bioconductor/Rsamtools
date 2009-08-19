#include "snap.h"
#include "encode.h"
#include "IRanges_interface.h"

typedef struct _snap_elt_t *_SNAP_ELT_PTR;

typedef struct _snap_elt_t {
	void *start, *curr, *end;
	_SNAP_ELT_PTR next;
} _SNAP_ELT_T;

typedef struct _snap_t {
	_SNAP_ELT_PTR char_root, char_curr;
	_SNAP_ELT_PTR width_root, width_curr;
} _SNAP_T;

const int _SNAP_BUF_SZ = 1048576;

_SNAP_ELT_PTR
_snap_elt_new(size_t buf_sz, int sizeof_elt)
{
	_SNAP_ELT_PTR selt = (_SNAP_ELT_PTR) R_alloc(1, sizeof(_SNAP_ELT_T));
	selt->start = selt->curr = (void *) R_alloc(buf_sz, sizeof_elt);
	selt->end = selt->start + buf_sz * sizeof_elt;
	selt->next = NULL;
	return selt;
}

_SNAP_PTR
_snap_new()
{
	_SNAP_PTR sptr = (_SNAP_PTR) R_alloc(1, sizeof(_SNAP_T));
	sptr->char_root = sptr->char_curr = 
		_snap_elt_new(_SNAP_BUF_SZ, sizeof(char));
	sptr->width_root = sptr->width_curr =
		_snap_elt_new(_SNAP_BUF_SZ, sizeof(int));
	return sptr;
}

_SNAP_ELT_PTR
_snap_add_elt(_SNAP_ELT_PTR ptr, int len, int sizeof_elt)
{
	_SNAP_ELT_PTR elt = _snap_elt_new(len, sizeof_elt);
	ptr->next = elt;
	return elt;
}

void
_snap_append(_SNAP_PTR sptr, const char *string)
{
	int len = strlen(string);
	if (sptr->char_curr->curr + len > sptr->char_curr->end) {
		size_t blen = len > _SNAP_BUF_SZ ? len : _SNAP_BUF_SZ;
		sptr->char_curr = 
			_snap_add_elt(sptr->char_curr, blen, sizeof(char));
	}
	memcpy(sptr->char_curr->curr, string, len);
	sptr->char_curr->curr += len;
	if (sptr->width_curr->curr == sptr->width_curr->end)
		sptr->width_curr =
			_snap_add_elt(sptr->width_curr, _SNAP_BUF_SZ, sizeof(int));
	*((int *) sptr->width_curr->curr) = len;
	sptr->width_curr->curr += sizeof(int);
}

SEXP
_snap_to_XStringSet(SEXP seq, SEXP width, const char *baseclass)
{
	ENCODE_FUNC encode = _encoder(baseclass);
	char *str = (char *) RAW(seq);
	for (int i = 0; i < LENGTH(seq); ++i)		
		str[i] = encode(str[i]);
	SEXP ptr = PROTECT(new_SequencePtr("RawPtr", seq));
	SEXP xstring = PROTECT(new_XSequence(baseclass, ptr, 0, LENGTH(ptr)));

	SEXP start = PROTECT(NEW_INTEGER(LENGTH(width)));
	int s = 1;
	for (int i = 0; i < LENGTH(width); ++i) {
		INTEGER(start)[i] = s;
		s += INTEGER(width)[i];
	}
	SEXP irange = 
		PROTECT(new_IRanges("IRanges", start, width, R_NilValue));

	const int XSETCLASS_BUF = 40;
	if (strlen(baseclass) > XSETCLASS_BUF - 4)
		error("Rsamtools internal error: *Set buffer too small");
	char xsetclass[XSETCLASS_BUF];
	snprintf(xsetclass, XSETCLASS_BUF, "%sSet", baseclass);
	SEXP xclassdef = PROTECT(MAKE_CLASS(xsetclass));
	SEXP xstringset = PROTECT(NEW_OBJECT(xclassdef));
	SET_SLOT(xstringset, PROTECT(mkChar("super")), xstring);
	SET_SLOT(xstringset, PROTECT(mkChar("ranges")), irange);

	UNPROTECT(8);
	return xstringset;
}

size_t
_snap_elts_length(_SNAP_ELT_PTR sptr)
{
	size_t len = 0;
	while (sptr != NULL) {
		len += sptr->curr - sptr->start;
		sptr = sptr->next;
	}
	return len;
}

void
_snap_elts_to(_SNAP_ELT_PTR sptr, void *to)
{
	while (sptr != NULL) {
		size_t len = sptr->curr - sptr->start;
		memcpy(to, sptr->start, len);
		to += len;
		sptr = sptr->next;
	}
}

SEXP
_snap_as_XStringSet(_SNAP_PTR sptr, const char *baseclass)
{
	size_t 
		char_len = _snap_elts_length(sptr->char_root) / sizeof(char),
		width_len = _snap_elts_length(sptr->width_root) / sizeof(int);
	SEXP
		r_char = PROTECT(NEW_RAW(char_len)),
		r_width = PROTECT(NEW_INTEGER(width_len));
	_snap_elts_to(sptr->char_root, RAW(r_char));
	_snap_elts_to(sptr->width_root, INTEGER(r_width));

	SEXP xstringset = _snap_to_XStringSet(r_char, r_width, baseclass);
	UNPROTECT(2);
	return xstringset;
}
