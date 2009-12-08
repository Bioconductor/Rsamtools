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
		SEXP str = STRING_ELT(cigar, i);
		if (NA_STRING == str) {
			n = 0;
		} else {
			const char *c = CHAR(str);
			while (*c != '\0') {
				if (57 < *c) n += 1;
				c += 1;
			}
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
		if (NA_STRING == str) {
			SET_STRING_ELT(cidx, j, str);
			INTEGER(element)[j] = NA_INTEGER;
			INTEGER(length)[j] = NA_INTEGER;
			SET_STRING_ELT(value, j, NA_STRING);
			j += 1;
		} else {
			const char *c = CHAR(str);
			int n = 0, elt = 1;
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
 * valid_cigar()
 * split_cigar()
 * cigar_to_qwidth()
 * cigar_to_IRanges()
 * cigar_to_GappedRanges()
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

static const char *cigar_string_to_qwidth(SEXP cigar_string, int clip_reads,
		int *qwidth)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	cig0 = CHAR(cigar_string);
	*qwidth = offset = 0;
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
		if (n == -1)
			return errmsg_buf;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': *qwidth += OPL; break;
		/* Insertion to the reference */
		    case 'I': *qwidth += OPL; break;
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N': break;
		/* Soft clip on the read */
		    case 'S': *qwidth += OPL; break;
		/* Hard clip on the read */
		    case 'H': if (!clip_reads) *qwidth += OPL; break;
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
		offset += n;
	}
	return NULL;
}

static const char *split_cigar_string(SEXP cigar_string,
		CharAE *OPbuf, IntAE *OPLbuf)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	cig0 = CHAR(cigar_string);
	offset = 0;
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
		if (n == -1)
			return errmsg_buf;
		CharAE_insert_at(OPbuf, OPbuf->nelt, OP);
		IntAE_insert_at(OPLbuf, OPLbuf->nelt, OPL);
		offset += n;
	}
	return NULL;
}

static void append_range(RangeAE *range_ae, int start, int width)
{
	RangeAE_insert_at(range_ae, range_ae->start.nelt, start, width);
}

static const char *cigar_string_to_ranges(SEXP cigar_string, int pos_elt,
		int Ds_as_Ns, RangeAE *range_ae)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */, start, width;
	char OP /* Operation */;

	cig0 = CHAR(cigar_string);
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
		/* Deletion from the reference */
		    case 'D':
			if (Ds_as_Ns)
				start += OPL;
			else
				width = OPL;
		    break;
		/* Skipped region from the reference */
		    case 'N': start += OPL; break;
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
		if (width) {
			append_range(range_ae, start, width);
			start += width;
		}
		offset += n;
	}
	return NULL;
}

/* Unlike cigar_string_to_ranges(), cigar_string_to_ranges2() merges adjacent
   ranges. */
static const char *cigar_string_to_ranges2(SEXP cigar_string, int pos_elt,
		int Ds_as_Ns, RangeAE *range_ae)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */, start, width;
	char OP /* Operation */;

	cig0 = CHAR(cigar_string);
	offset = 0;
	start = pos_elt;
	width = 0;
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
		if (n == -1)
			return errmsg_buf;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': width += OPL; break;
		/* Insertion to the reference */
		    case 'I': break;
		/* Deletion from the reference */
		    case 'D':
			if (Ds_as_Ns) {
				if (width)
					append_range(range_ae, start, width);
				start += width + OPL;
				width = 0;
			} else {
				width += OPL;
			}
		    break;
		/* Skipped region from the reference */
		    case 'N':
			if (width)
				append_range(range_ae, start, width);
			start += width + OPL;
			width = 0;
		    break;
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
		offset += n;
	}
	if (width)
		append_range(range_ae, start, width);
	return NULL;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 *   ans_type: a single integer specifying the type of answer to return:
 *     0: 'ans' is a string describing the first validity failure or NULL;
 *     1: 'ans' is logical vector with TRUE values for valid elements
 *        in 'cigar'.
 */
SEXP valid_cigar(SEXP cigar, SEXP ans_type)
{
	SEXP ans, cigar_string;
	int cigar_length, ans_type0, i, qwidth;
	const char *errmsg;
	char string_buf[200];

	cigar_length = LENGTH(cigar);
	ans_type0 = INTEGER(ans_type)[0];
	if (ans_type0 == 1)
		PROTECT(ans = NEW_LOGICAL(cigar_length));
	else
		ans = R_NilValue;
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		/* we use cigar_string_to_qwidth() here just for its ability
                   to parse and detect ill-formed CIGAR strings */
		errmsg = cigar_string_to_qwidth(cigar_string, 1, &qwidth);
		if (ans_type0 == 1) {
			LOGICAL(ans)[i] = errmsg == NULL;
			continue;
		}
		if (errmsg != NULL) {
			snprintf(string_buf, sizeof(string_buf),
				 "element %d is invalid (%s)", i + 1, errmsg);
			return mkString(string_buf);
		}
	}
	if (ans_type0 == 1)
		UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read.
 * Return a list of the same length as 'cigar' where each element is itself
 * a list with 2 elements of the same lengths, the 1st one being a raw
 * vector containing the CIGAR operations and the 2nd one being an integer
 * vector containing the lengths of the CIGAR operations.
 */
SEXP split_cigar(SEXP cigar)
{
	SEXP ans, cigar_string, ans_elt, ans_elt_elt0, ans_elt_elt1;
	int cigar_length, i;
	CharAE OPbuf;
	IntAE OPLbuf;
	const char *errmsg;

	cigar_length = LENGTH(cigar);
	PROTECT(ans = NEW_LIST(cigar_length));
	OPbuf = new_CharAE(0);
	OPLbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			UNPROTECT(1);
			error("'cigar' contains NAs");
		}
		OPbuf.nelt = OPLbuf.nelt = 0;
		errmsg = split_cigar_string(cigar_string, &OPbuf, &OPLbuf);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		PROTECT(ans_elt = NEW_LIST(2));
		PROTECT(ans_elt_elt0 = CharAE_asRAW(&OPbuf));
		PROTECT(ans_elt_elt1 = IntAE_asINTEGER(&OPLbuf));
		SET_VECTOR_ELT(ans_elt, 0, ans_elt_elt0);
		SET_VECTOR_ELT(ans_elt, 1, ans_elt_elt1);
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(3);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 *   before_hard_clipping: TRUE or FALSE indicating whether the returned
 *          widths should be those of the reads before or after "hard
 *          clipping".
 * Return an integer vector of the same length as 'cigar' containing the
 * widths of the reads as inferred from the cigar information.
 */
SEXP cigar_to_qwidth(SEXP cigar, SEXP before_hard_clipping)
{
	SEXP ans, cigar_string;
	int clip_reads, cigar_length, i, qwidth;
	const char *errmsg;

	clip_reads = !LOGICAL(before_hard_clipping)[0];
	cigar_length = LENGTH(cigar);
	PROTECT(ans = NEW_INTEGER(cigar_length));
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			INTEGER(ans)[i] = NA_INTEGER;
			continue;
		}
		errmsg = cigar_string_to_qwidth(cigar_string, clip_reads,
				&qwidth);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		INTEGER(ans)[i] = qwidth;
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character string containing the extended CIGAR;
 *   drop_D_ranges: TRUE or FALSE indicating whether Ds should be treated
 *          like Ns or not;
 *   merge_ranges: TRUE or FALSE indicating whether adjacent ranges coming
 *          from the same cigar should be merged or not.
 * Return an IRanges object describing the alignment.
 */
SEXP cigar_to_IRanges(SEXP cigar, SEXP drop_D_ranges, SEXP merge_ranges)
{
	RangeAE range_ae;
	SEXP cigar_string;
	int Ds_as_Ns, merge_ranges0;
	const char *errmsg;

	cigar_string = STRING_ELT(cigar, 0);
	if (cigar_string == NA_STRING)
		error("'cigar' is NA");
	Ds_as_Ns = LOGICAL(drop_D_ranges)[0];
	merge_ranges0 = LOGICAL(merge_ranges)[0];
	range_ae = new_RangeAE(0, 0);
	errmsg = merge_ranges0 ?
			cigar_string_to_ranges2(cigar_string, 1,
				Ds_as_Ns, &range_ae) :
			cigar_string_to_ranges(cigar_string, 1,
				Ds_as_Ns, &range_ae);
	if (errmsg != NULL)
		error("%s", errmsg);
	return RangeAE_asIRanges(&range_ae);
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 *   pos:   integer vector containing the 1-based leftmost position/coordinate
 *          of the clipped read sequence;
 *   flag:  NULL or an integer vector containing the SAM flag for each
 *          read;
 *   drop_D_ranges: TRUE or FALSE indicating whether Ds should be treated
 *          like Ns or not;
 * 'cigar', 'pos' and 'flag' (when not NULL) are assumed to have the same
 * length (which is the number of aligned reads).
 *
 * Returns a GappedRanges object of the same length as the input.
 * NOTE: See note for cigar_to_list_of_IRanges() below about the strand.
 * TODO: Support character factor 'cigar' in addition to current character
 *       vector format.
 */
SEXP cigar_to_GappedRanges(SEXP cigar, SEXP pos, SEXP flag, SEXP drop_D_ranges)
{
	SEXP cigar_string, ans, ans_cirl, ans_cirl_unlistData,
	     ans_cirl_partitioning, ans_cirl_partitioning_end;
	int cigar_length, Ds_as_Ns, i, pos_elt, flag_elt;
	RangeAE range_ae;
	const char *errmsg;

	cigar_length = LENGTH(cigar);
	Ds_as_Ns = LOGICAL(drop_D_ranges)[0];
	/* we will generate at least 'cigar_length' ranges, and possibly more */
	range_ae = new_RangeAE(cigar_length, 0);
	PROTECT(ans_cirl_partitioning_end = NEW_INTEGER(cigar_length));
	for (i = 0; i < cigar_length; i++) {
		if (flag != R_NilValue) {
			flag_elt = INTEGER(flag)[i];
			if (flag_elt == NA_INTEGER) {
				UNPROTECT(1);
				error("'flag' contains NAs");
			}
			if (flag_elt & 0x404)
				continue;
		}
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			UNPROTECT(1);
			error("'cigar' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		}
		pos_elt = INTEGER(pos)[i];
		if (pos_elt == NA_INTEGER) {
			UNPROTECT(1);
			error("'pos' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		}
		errmsg = cigar_string_to_ranges2(cigar_string, pos_elt,
				Ds_as_Ns, &range_ae);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		INTEGER(ans_cirl_partitioning_end)[i] = range_ae.start.nelt;
	}
	// TODO: Add C-level constructors in IRanges for PartitioningByEnd,
	// CompressedIRangesList and GappedRanges objects and use them here.
	PROTECT(ans_cirl_unlistData = RangeAE_asIRanges(&range_ae));
	PROTECT(ans_cirl_partitioning = NEW_OBJECT(MAKE_CLASS("PartitioningByEnd")));
	SET_SLOT(ans_cirl_partitioning, install("end"), ans_cirl_partitioning_end);
	PROTECT(ans_cirl = NEW_OBJECT(MAKE_CLASS("CompressedIRangesList")));
	SET_SLOT(ans_cirl, install("unlistData"), ans_cirl_unlistData);
	SET_SLOT(ans_cirl, install("partitioning"), ans_cirl_partitioning);
	PROTECT(ans = NEW_OBJECT(MAKE_CLASS("GappedRanges")));
	SET_SLOT(ans, install("cirl"), ans_cirl);
	UNPROTECT(5);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 *   rname: character factor containing the name of the reference sequence
 *          associated with each read (i.e. the name of the sequence the
 *          read has been aligned to);
 *   pos:   integer vector containing the 1-based leftmost position/coordinate
 *          of the clipped read sequence;
 *   flag:  NULL or an integer vector containing the SAM flag for each
 *          read;
 *   drop_D_ranges: TRUE or FALSE indicating whether Ds should be treated
 *          like Ns or not;
 *   merge_ranges: TRUE or FALSE indicating whether adjacent ranges coming
 *          from the same cigar should be merged or not.
 * 'cigar', 'rname', 'pos' and 'flag' (when not NULL) are assumed to have
 * the same length (which is the number of aligned reads).
 *
 * Return a list of IRanges objects named with the factor levels in 'rname'.
 *
 * NOTE: According to the SAM Format Specification (0.1.2-draft 20090820),
 *   the CIGAR (and the read sequence) stored in the SAM file are represented
 *   on the + strand of the reference sequence. This means that, for a read
 *   that aligns to the - strand, the bases have been reverse complemented
 *   from the unmapped read sequence, and that the corresponding CIGAR has
 *   been reversed. So it seems that, for now, we don't need to deal with the
 *   strand information at all (as long as we are only interested in
 *   returning a list of IRanges objects that is suitable for coverage
 *   extraction).
 *
 * TODO:
 * - Support 'rname' of length 1.
 * - Support character factor 'cigar' in addition to current character vector
 *   format.
 */
SEXP cigar_to_list_of_IRanges(SEXP cigar, SEXP rname, SEXP pos,
		SEXP flag, SEXP drop_D_ranges, SEXP merge_ranges)
{
	SEXP rname_levels, cigar_string, ans, ans_names;
	int ans_length, nreads, Ds_as_Ns, merge_ranges0,
	    i, level, pos_elt, flag_elt;
	RangeAEAE range_aeae;
	const char *errmsg;

	rname_levels = GET_LEVELS(rname);
	ans_length = LENGTH(rname_levels);
	range_aeae = new_RangeAEAE(ans_length, ans_length);
	nreads = LENGTH(pos);
	Ds_as_Ns = LOGICAL(drop_D_ranges)[0];
	merge_ranges0 = LOGICAL(merge_ranges)[0];
	for (i = 0; i < nreads; i++) {
		if (flag != R_NilValue) {
			flag_elt = INTEGER(flag)[i];
			if (flag_elt == NA_INTEGER)
				error("'flag' contains NAs");
			if (flag_elt & 0x404)
				continue;
		}
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING)
			error("'cigar' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		level = INTEGER(rname)[i];
		if (level == NA_INTEGER)
			error("'rname' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		pos_elt = INTEGER(pos)[i];
		if (pos_elt == NA_INTEGER)
			error("'pos' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		errmsg = merge_ranges0 ?
			cigar_string_to_ranges2(cigar_string, pos_elt,
				Ds_as_Ns, range_aeae.elts + level - 1) :
			cigar_string_to_ranges(cigar_string, pos_elt,
				Ds_as_Ns, range_aeae.elts + level - 1);
		if (errmsg != NULL)
			error("in 'cigar' element %d: %s", i + 1, errmsg);
	}
	PROTECT(ans = RangeAEAE_asLIST(&range_aeae));
	PROTECT(ans_names = duplicate(rname_levels));
	SET_NAMES(ans, ans_names);
	UNPROTECT(2);
	return ans;
}

