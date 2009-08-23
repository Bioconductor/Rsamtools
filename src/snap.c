#include "snap.h"
#include "encode.h"
#include "utilities.h"
#include "IRanges_interface.h"

/* 
   An attempt at better memory management. For strings (e.g.,
   sequences or base qualities), allocate relatively large blocks of
   memory in a linked list. As the data grows, memory gets consumed
   without copying. The final stage transforms the stored data into R
   structures, requiring one complete copy of the string buffers and
   of the corresponding widths.
*/

/* _SNAP_ELT_T: a block of memory */

const int _SNAP_ELT_SZ = 1048576;

typedef struct _snap_elt_t _SNAP_ELT_T;

struct _snap_elt_t {
	void *start, *curr, *end;
	_SNAP_ELT_T *next;
};

#define SE_ADDR_BASE(snapelt_p) ((snapelt_p)->start)
#define SE_BYTES_AVAIL(snapelt_p) \
	((snapelt_p)->end - (snapelt_p)->curr)
#define SE_BYTES_USED(snapelt_p) \
	((snapelt_p)->curr - SE_ADDR_BASE(snapelt_p))

/* _SNAP_LIST_T: a dynamic linked list of memory blocks */

typedef struct _snap_list_t _SNAP_LIST_T;

struct _snap_list_t {
	SNAPTYPE type;
	size_t sizeof_type;
	_SNAP_ELT_T *root, *curr;
};

size_t	_snap_list_byte_length(_SNAP_LIST_T *);
void _snap_list_delete(_SNAP_LIST_T *);

#define SL_TYPE(snaplist_p) ((snaplist_p)->type)
#define SL_SIZEOF_TYPE(snaplist_p) ((snaplist_p)->sizeof_type)
#define SL_ADDR(snaplist_p) ((snaplist_p)->curr->curr)
#define SL_LENGTH(snaplist_p)  \
	(_snap_list_byte_length(snaplist_p) / SL_SIZEOF_TYPE(snaplist_p))
#define SL_BYTES_AVAIL(snaplist_p) \
	(SE_BYTES_AVAIL((snaplist_p)->curr))

/* _SNAP_T: a collection of memory blocks, one for character data, one
 * for integer */

struct _snap_t {
	_SNAP_LIST_T *str, *width;
};

/* register snaps here, so they can be recovered in case of error */

struct _snap_alloc_t {
	_SNAP_T *snap;
	struct _snap_alloc_t *next;
};

static struct _snap_alloc_t *_snap_alloc_g = NULL;

void
_snap_alloc_register(_SNAP_T *snap)
{
	struct _snap_alloc_t *node = Calloc(1, struct _snap_alloc_t);
	node->snap = snap;
	node->next = _snap_alloc_g;
	_snap_alloc_g = node;
}

void
__snap_free(_SNAP_T *snap)
{
	_snap_list_delete(snap->str);
	_snap_list_delete(snap->width);
	Free(snap);
}

void
_snap_delete(_SNAP_T *snap)
{
	struct _snap_alloc_t *prev = NULL, *node = _snap_alloc_g;
	while (node != NULL) {
		if (node->snap == snap) {
			if (prev == NULL)
				_snap_alloc_g = node->next;
			else
				prev->next = node->next;
			__snap_free(node->snap);
			Free(node);
			break;
		}
		prev = node;
		node = node->next;
	}
}

/* only called during error handling */
void
_snap_cleanup_all()
{
	struct _snap_alloc_t *node = _snap_alloc_g;
	while (node != NULL) {
		__snap_free(node->snap);
		struct _snap_alloc_t *this = node;
		node = node->next;
		Free(this);
	}		
	_snap_alloc_g = NULL;
}

/* implementation */

_SNAP_ELT_T *
_snap_elt_new(SNAPTYPE type, size_t elt_len)
{
	_SNAP_ELT_T *selt = Calloc(1, _SNAP_ELT_T);
	void *buf = 0;
	size_t sizeof_type = 0;
	switch(type) {
	case INTTYPE:
		buf = (void *) Calloc(elt_len, int);
		sizeof_type = sizeof(int);
		break;
	case RAWTYPE:
		buf = (void *) Calloc(elt_len, char);
		sizeof_type = sizeof(char);
		break;
	default:
		Rf_error("Rsamtools invalid _snap_list_new type '%d'",
				 type);
	}
	selt->start = selt->curr = buf;
	selt->end = selt->start + elt_len * sizeof_type;
	selt->next = NULL;
	return selt;
}

_SNAP_LIST_T *
_snap_list_new(SNAPTYPE type, size_t sizeof_elt)
{
	_SNAP_LIST_T *slptr = Calloc(1, _SNAP_LIST_T);
	slptr->type = type;
	switch(type) {
	case INTTYPE:
		slptr->sizeof_type = sizeof(int);
		break;
	case RAWTYPE:
		slptr->sizeof_type = sizeof(char);
		break;
	default:
		Rf_error("Rsamtools invalid _snap_list_new type '%d'",
				 type);
	}
	slptr->root = slptr->curr = 
		_snap_elt_new(type, sizeof_elt);
	return slptr;
}

void
_snap_list_delete(_SNAP_LIST_T *lst)
{
	_SNAP_ELT_T *elt = lst->root;
	while (elt != NULL) {
		_SNAP_ELT_T *tmp = elt->next;
		Free(elt);
		elt = tmp;
	}
	Free(lst);
}

void
_snap_list_add_elt(_SNAP_LIST_T *ptr, size_t len)
{
	_SNAP_ELT_T *elt = 
		_snap_elt_new(SL_TYPE(ptr), len);
	ptr->curr->next = elt;
	ptr->curr = elt;
}

size_t
_snap_list_byte_length(_SNAP_LIST_T *lst)
{
	_SNAP_ELT_T *ptr = lst->root;
	size_t len = 0;
	while (ptr != NULL) {
		len += SE_BYTES_USED(ptr);
		ptr = ptr->next;
	}
	return len;
}

void
_snap_list_append_v(_SNAP_LIST_T *lst, const void *data, R_len_t n_data_elts)
{
	/* following use bytes, not elements */
	size_t len = SL_SIZEOF_TYPE(lst) * n_data_elts;
	if (SL_BYTES_AVAIL(lst) < len) {
		size_t buflen = len > _SNAP_ELT_SZ ? len : _SNAP_ELT_SZ;
		_snap_list_add_elt(lst, buflen);
	}
	memcpy(SL_ADDR(lst), data, len);
	SL_ADDR(lst) += len;
}

void
_snap_list_append_int(_SNAP_LIST_T *lst, int x)
{
	size_t len = SL_SIZEOF_TYPE(lst);
	if (SL_BYTES_AVAIL(lst) < len)
		_snap_list_add_elt(lst, _SNAP_ELT_SZ);
	*((int *) SL_ADDR(lst)) = x;
	SL_ADDR(lst) += len;
}

/* snap: higher level interface */

_SNAP_T *
_snap_new()
{
	_SNAP_T *sptr = Calloc(1, _SNAP_T);
	sptr->str = _snap_list_new(RAWTYPE, _SNAP_ELT_SZ);
	sptr->width = _snap_list_new(INTTYPE, _SNAP_ELT_SZ);
	_snap_alloc_register(sptr);
	return sptr;
}

void
_snap_append(_SNAP_T *sptr, const char *string)
{
	int len = strlen(string);
	_snap_list_append_v(sptr->str, string, len);
	_snap_list_append_int(sptr->width, len);
}

/* collapse into R objects */

SEXP
_snap_list_to(_SNAP_LIST_T *lst)
{
	SNAPTYPE type = SL_TYPE(lst);
	size_t n_elt = SL_LENGTH(lst);
	SEXP sexp = 0;
	void *to = 0;
	switch(type) {
	case RAWTYPE:
		sexp = PROTECT(NEW_RAW(n_elt));
		to = RAW(sexp);
		break;
	case INTTYPE:
		sexp = PROTECT(NEW_INTEGER(n_elt));
		to = INTEGER(sexp);
		break;
	default:
		Rf_error("Rsamtools invalid _snap_list_to type '%d'",
				 type);
	}

	_SNAP_ELT_T *ptr = lst->root;
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
	SEXP xstring = PROTECT(new_XSequence(baseclass, ptr, 0, LENGTH(seq)));

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
        
	SET_SLOT(xstringset, install("super"), xstring);
	SET_SLOT(xstringset, install("ranges"), irange);

	UNPROTECT(6);
	return xstringset;
}

SEXP
_snap_as_PhredQuality(_SNAP_T *snap)
{
	SEXP xstringset = PROTECT(_snap_as_XStringSet(snap, "BString"));
	SEXP s, t, nmspc, result;
	nmspc = PROTECT(_get_namespace("Rsamtools"));
	NEW_CALL(s, t, "PhredQuality", nmspc, 2);
	CSET_CDR(t, "x", xstringset);
	CEVAL_TO(s, nmspc, result);
	UNPROTECT(2);
	return result;
}

SEXP
_snap_as_XStringSet(_SNAP_T *sptr, const char *baseclass)
{
	SEXP r_str = PROTECT(_snap_list_to(sptr->str)),
		r_width = PROTECT(_snap_list_to(sptr->width));
	SEXP xstringset = _snap_to_XStringSet(r_str, r_width, baseclass);
	UNPROTECT(2);
	return xstringset;
}
