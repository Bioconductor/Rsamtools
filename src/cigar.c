#include <Rdefines.h>
#include "cigar.h"
#include "utilities.h"
#include "IRanges_interface.h"

#include <ctype.h> /* for isdigit() */

SEXP
cigar_run_count(SEXP cigar)
{
	SEXP count = NEW_INTEGER(LENGTH(cigar));
	for (int i = 0; i < LENGTH(cigar); ++i) {
		int n = 0;
		const char *c = CHAR(STRING_ELT(cigar, i));
		while (*c != '\0') {
			if (57 < *c) n += 1;
			c += 1;
		}
		INTEGER(count)[i] = n;
	}
	return count;
}

SEXP
cigar_table(SEXP cigar)
{
	int i, j, len = 0;
	SEXP rcount = PROTECT(cigar_run_count(cigar));
	for (i = 0; i < LENGTH(rcount); ++i)
		if (INTEGER(rcount)[i] == 0)
			len += 1; /* single-line entry for zero-length cigar "" */
		else 
			len += INTEGER(rcount)[i];
	SEXP cidx = PROTECT(NEW_CHARACTER(len)),
		element = PROTECT(NEW_INTEGER(len)),
		length = PROTECT(NEW_INTEGER(len)),
		value = PROTECT(NEW_CHARACTER(len));

	for (i = j = 0; i < LENGTH(cigar); ++i) {
		SEXP str = STRING_ELT(cigar, i);
		const char *c = CHAR(str);
		int n = 0, elt = 1;
		if (*c == '\0') {		/* zero-length cigar */
			SET_STRING_ELT(cidx, j, str);
			INTEGER(element)[j] = NA_INTEGER;
			INTEGER(length)[j] = NA_INTEGER;
			SET_STRING_ELT(value, j, NA_STRING);
			j += 1;
		} else {
			while (*c != '\0') {
				if (57 >= *c) {
					n = n * 10 + (int) (*c - 48);
				} else {
					SET_STRING_ELT(cidx, j, str);
					INTEGER(element)[j] = elt;
					INTEGER(length)[j] = n;
					SET_STRING_ELT(value, j, mkCharLen(c, 1));
					n = 0; 
					elt += 1;
					j += 1;
				}
				c += 1;
			}
		}
	}

	/* construct the table */
	SEXP s, t, nmspc, result;
    SEXP m_false = PROTECT(allocVector(LGLSXP, 1));
    LOGICAL(m_false)[0] = 0;
	nmspc = PROTECT(_get_namespace("Rsamtools"));
	NEW_CALL(s, t, "data.frame", nmspc, 6);
	CSET_CDR(t, "cigar", cidx);
	CSET_CDR(t, "element", element);
	CSET_CDR(t, "length", length);
	CSET_CDR(t, "value", value);
	CSET_CDR(t, "stringsAsFactors", m_false);
	CEVAL_TO(s, nmspc, result);
	UNPROTECT(7);
	return result;
}


/****************************************************************************
 * cigar_to_list_of_IRanges()
 */

static char errmsg_buf[200];

/* Return the number of chars that was read, or 0 if there is no more char
   to read (i.e. cig0[offset] is '\0'), or -1 in case of a parse error. */
/*
static int get_next_cigar_OP(const char *cig0, int offset,
		int *OPL, char *OP)
{
	char c;
	int ret, n;

	if (!(c = cig0[offset]))
		return 0;
	if (!isdigit(c)) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "unsigned decimal integer expected at char %d",
			 offset + 1);
		return -1;
	}
	ret = sscanf(cig0 + offset, "%d%c%n", OPL, OP, &n);
	if (ret < 2) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "CIGAR parse error at or after char %d",
			 offset + 1);
		return -1;
	}
	if (*OPL <= 0) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "invalid CIGAR operation length at char %d",
			 offset + 1);
		return -1;
	}
	return n;
}
*/

static int get_next_cigar_OP(const char *cig0, int offset,
		int *OPL, char *OP)
{
	char c;
	int offset0, opl;

	if (!(c = cig0[offset]))
		return 0;
	offset0 = offset;
	opl = 0;
	while (isdigit(c = cig0[offset])) {
		offset++;
		opl *= 10;
		opl += c - '0';
	}
	if (opl == 0) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "invalid CIGAR operation length at char %d",
			 offset + 1);
		return -1;
	}
	*OPL = opl;
	if (!(*OP = cig0[offset])) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "unexpected CIGAR end at char %d",
			 offset + 1);
		return -1;
	}
	offset++;
	return offset - offset0;
}

static const char *expand_cigar(RangeAE *range_ae, int pos_elt, SEXP cigar_elt)
{
	const char *cig0;
	int offset, OPL /* Operation Length */, n, start, width;
	char OP /* Operation */;

	cig0 = CHAR(cigar_elt);
	offset = 0;
	start = pos_elt;
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
		if (n == -1)
			return errmsg_buf;
		width = 0;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': width = OPL; break;
		/* Insertion to the reference */
		    case 'I': break;
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N': start += OPL; break;
		/* Soft/hard clip on the read */
		    case 'S': case 'H': break;
		/* Padding (silent deletion from the padded reference
		   sequence) */
		    case 'P':
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "CIGAR operation '%c' (at char %d) is not "
				 "supported yet, sorry!", OP, offset + 1);
			return errmsg_buf;
		    break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		if (width != 0) {
			RangeAE_insert_at(range_ae, range_ae->start.nelt,
					  start, width);
			start += width;
		}
		offset += n;
	}
	return NULL;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   rname: character factor containing the name of the reference sequence
 *          associated with each read (i.e. the name of the sequence the
 *          read has been aligned to);
 *   strand: ignored for now;
 *   pos: integer vector containing the 1-based leftmost position/coordinate
 *          of the clipped read sequence;
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 * 'rname', 'pos' and 'cigar' are assumed to have the same length (which is
 * the number of aligned reads).
 *
 * Return a list of IRanges objects named with the factor levels in 'rname'.
 *
 * TODO:
 * - Decide what we want to do with the 'strand' arg. Our current
 *   understanding of the SAM Format Specification (0.1.2-draft 20090820)
 *   is that the CIGAR (and the read sequence) stored in the SAM file are
 *   represented on the + strand of the reference sequence. This means that,
 *   for a read that aligns to the - strand, the bases have been reverse
 *   complemented from the unmapped read sequence, and that the corresponding
 *   CIGAR has been reversed.
 *   So it seems that, for now, we don't need to deal with the strand
 *   information at all (as long as we are only interested in returning a list
 *   of IRanges objects that is suitable for coverage extraction).
 * - Support 'rname' of length 1.
 * - Support character factor 'cigar' in addition to current character vector
 *   format.
 */
SEXP cigar_to_list_of_IRanges(SEXP rname, SEXP strand, SEXP pos, SEXP cigar)
{
	SEXP rname_levels, cigar_elt, ans, ans_names;
	int ans_length, nreads, i, level, pos_elt;
	RangeAEAE range_aeae;
	const char *errmsg;

	rname_levels = GET_LEVELS(rname);
	ans_length = LENGTH(rname_levels);
	range_aeae = new_RangeAEAE(ans_length, ans_length);
	nreads = LENGTH(pos);
	for (i = 0; i < nreads; i++) {
		level = INTEGER(rname)[i];
		if (level == NA_INTEGER)
			error("'rname' contains NAs");
		pos_elt = INTEGER(pos)[i];
		if (pos_elt == NA_INTEGER)
			error("'pos' contains NAs");
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING)
			error("'cigar' contains NAs");
		errmsg = expand_cigar(range_aeae.elts + level - 1,
				      pos_elt, cigar_elt);
		if (errmsg != NULL)
			error("in 'cigar' element %d: %s", i + 1, errmsg);
	}
	PROTECT(ans = RangeAEAE_asLIST(&range_aeae));
	PROTECT(ans_names = duplicate(rname_levels));
	SET_NAMES(ans, ans_names);
	UNPROTECT(2);
	return ans;
}

