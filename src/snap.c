#include "snap.h"
#include "encode.h"
#include "IRanges_interface.h"

/* _SNAP_ELT_T: a block of memory to be filled */

const int _SNAP_ELT_SZ = 1048576;

typedef struct _snap_elt_t *_SNAP_ELT_PTR;
typedef struct _snap_elt_t {
	void *start, *curr, *end;
	_SNAP_ELT_PTR next;
} _SNAP_ELT_T;

#define SE_ADDR_BASE(snapelt_p) ((snapelt_p)->start)
#define SE_BYTES_AVAIL(snapelt_p) \
	((snapelt_p)->end - (snapelt_p)->curr)
#define SE_BYTES_USED(snapelt_p) \
	((snapelt_p)->curr - SE_ADDR_BASE(snapelt_p))

/* _SNAP_LIST_T: a dynamic linked list of memory blocks */

typedef struct _snap_list_t *_SNAP_LIST_PTR;
typedef struct _snap_list_t {
	SEXPTYPE stype;
	_SNAP_ELT_PTR root, curr;
} _SNAP_LIST_T;

size_t _snap_list_elt_size(SEXPTYPE);
size_t	_snap_list_byte_length(_SNAP_LIST_PTR);

#define SL_TYPE(snaplist_p) ((snaplist_p)->stype)
#define SL_ADDR(snaplist_p) ((snaplist_p)->curr->curr)
#define SL_BYTE_LENGTH(snaplist_p) \
	(_snap_list_byte_length(snaplist_p))
#define SL_LENGTH(snaplist_p)  \
	(SL_BYTE_LENGTH(snaplist_p) / SL_ELT_SIZE(snaplist_p))
#define SL_ELT_SIZE(snaplist_p) \
	(_snap_list_elt_size(SL_TYPE(snaplist_p)))
#define SL_BYTES_AVAIL(snaplist_p) \
	(SE_BYTES_AVAIL((snaplist_p)->curr))

/* _SNAP_T: a collection of memory blocks, one for character data, one
 * for integer */
typedef struct _snap_t {
	_SNAP_LIST_T *str, *width;
} _SNAP_T;

_SNAP_ELT_PTR
_snap_elt_new(SEXPTYPE type, size_t sizeof_elt)
{
	int sizeof_i = _snap_list_elt_size(type);
	_SNAP_ELT_PTR selt = 
		(_SNAP_ELT_PTR) R_alloc(1, sizeof(_SNAP_ELT_T));
	selt->start = selt->curr = (void *) R_alloc(sizeof_elt, sizeof_i);
	selt->end = selt->start + sizeof_elt * sizeof_i;
	selt->next = NULL;
	return selt;
}

_SNAP_LIST_PTR
_snap_list_new(SEXPTYPE type, size_t sizeof_elt)
{
	_SNAP_LIST_PTR slptr = 
		(_SNAP_LIST_PTR) R_alloc(1, sizeof(_SNAP_LIST_T));
	slptr->stype = type;
	slptr->root = slptr->curr = _snap_elt_new(type, sizeof_elt);
	return slptr;
}

size_t
_snap_list_elt_size(SEXPTYPE type)
{
	size_t sz = 0;
	switch(type) {
	case INTSXP:
		sz = sizeof(int);
		break;
	case RAWSXP:
		sz = sizeof(char);
		break;
	default:
		Rf_error("Rsamtools invalid _snap_list_elt_size type '%d'",
				 type);
	}
	return sz;
}

void
_snap_list_add_elt(_SNAP_LIST_PTR ptr, size_t len)
{
	_SNAP_ELT_PTR elt = _snap_elt_new(SL_TYPE(ptr), len);
	ptr->curr->next = elt;
	ptr->curr = elt;
}

size_t
_snap_list_byte_length(_SNAP_LIST_PTR lst)
{
	_SNAP_ELT_PTR ptr = lst->root;
	size_t len = 0;
	while (ptr != NULL) {
		len += SE_BYTES_USED(ptr);
		ptr = ptr->next;
	}
	return len;
}

void
_snap_list_append_v(_SNAP_LIST_PTR lst, const void *data, R_len_t n_data_elts)
{
	/* following use bytes, not elements */
	size_t len = SL_ELT_SIZE(lst) * n_data_elts;
	if (SL_BYTES_AVAIL(lst) < len) {
		size_t buflen = len > _SNAP_ELT_SZ ? len : _SNAP_ELT_SZ;
		_snap_list_add_elt(lst, buflen);
	}
	memcpy(SL_ADDR(lst), data, len);
	SL_ADDR(lst) += len;
}

void
_snap_list_append_int(_SNAP_LIST_PTR lst, int x)
{
	if (SL_BYTES_AVAIL(lst) < SL_ELT_SIZE(lst))
		_snap_list_add_elt(lst, _SNAP_ELT_SZ);
	*((int *) SL_ADDR(lst)) = x;
	SL_ADDR(lst) += SL_ELT_SIZE(lst);
}

/* snap: higher level interface */

_SNAP_PTR
_snap_new()
{
	_SNAP_PTR sptr = (_SNAP_PTR) R_alloc(1, sizeof(_SNAP_T));
	sptr->str = _snap_list_new(RAWSXP, _SNAP_ELT_SZ);
	sptr->width = _snap_list_new(INTSXP, _SNAP_ELT_SZ);
	return sptr;
}

void
_snap_append(_SNAP_PTR sptr, const char *string)
{
	int len = strlen(string);
	_snap_list_append_v(sptr->str, string, len);
	_snap_list_append_int(sptr->width, len);
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

/* collapse into R objects */

SEXP
_snap_list_to(_SNAP_LIST_PTR lst)
{
	SEXPTYPE type = SL_TYPE(lst);
	size_t n_elt = SL_LENGTH(lst);
	SEXP sexp = PROTECT(allocVector(type, n_elt));
	void *to = 0;
	switch(type) {
	case RAWSXP:
		to = RAW(sexp);
		break;
	case INTSXP:
		to = INTEGER(sexp);
		break;
	default:
		Rf_error("Rsamtools invalid _snap_list_to type '%d'",
				 type);
	}

	_SNAP_ELT_PTR ptr = lst->root;
	while (ptr != NULL) {
		size_t len = SE_BYTES_USED(ptr);
		memcpy(to, SE_ADDR_BASE(ptr), len);
		to += len;
		ptr = ptr->next;
	}
	UNPROTECT(1);
	return sexp;
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

SEXP
_snap_as_XStringSet(_SNAP_PTR sptr, const char *baseclass)
{
	SEXP r_str = PROTECT(_snap_list_to(sptr->str)),
		r_width = PROTECT(_snap_list_to(sptr->width));
	SEXP xstringset = _snap_to_XStringSet(r_str, r_width, baseclass);
	UNPROTECT(2);
	return xstringset;
}
