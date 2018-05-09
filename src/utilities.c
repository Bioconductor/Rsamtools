#include <Rdefines.h>
#include <R_ext/RS.h>
#include "utilities.h"
#include "IRanges_interface.h"
#include "XVector_interface.h"
#include <htslib/khash.h>

void *_Rs_Realloc_impl(void *p, size_t n, size_t t)
{
    /* Realloc(p, 0, *) fails inappropriately */
    if (n == 0) {
        Free(p);
        p = NULL;
    } else {
        p = R_chk_realloc((void *) p, (size_t) (n * t));
    }
    return p;
}

SEXP _get_namespace(const char *pkg)
{
    SEXP fun = PROTECT(findFun(install("getNamespace"), R_GlobalEnv));
    SEXP nmspc = PROTECT(NEW_CHARACTER(1));
    SET_STRING_ELT(nmspc, 0, mkChar(pkg));
    nmspc = eval(lang2(fun, nmspc), R_GlobalEnv);
    UNPROTECT(2);
    return nmspc;
}

void _as_strand(SEXP vec)
{
    SEXP nmspc = PROTECT(_get_namespace("Rsamtools"));
    SEXP lvls = PROTECT(eval(findVar(install(".STRAND_LEVELS"), nmspc), nmspc));
    _as_factor_SEXP(vec, lvls);
    UNPROTECT(2);
}

void _as_nucleotide(SEXP vec)
{
    SEXP nmspc = PROTECT(_get_namespace("Rsamtools"));
    SEXP lvls = PROTECT(eval(findVar(install(".PILEUP_NUCLEOTIDE_LEVELS"),
                                     nmspc), nmspc));
    _as_factor_SEXP(vec, lvls);
    UNPROTECT(2);
}

SEXP _get_encoding_lookup(const char *from, const char *to)
{
    SEXP nmspc, fun, f, t, call, ans;
    nmspc = PROTECT(_get_namespace("Biostrings"));
    fun = findFun(install("get_seqtype_conversion_lookup"), nmspc);
    f = PROTECT(mkString(from));
    t = PROTECT(mkString(to));
    call = PROTECT(lang3(fun, f, t));
    ans = eval(call, nmspc);
    UNPROTECT(4);
    return ans;
}

SEXP _get_lkup(const char *baseclass)
{
    SEXP lkup = R_NilValue;
    switch (*baseclass) {
    case 'D':
        lkup = _get_encoding_lookup("B", "DNA");
        break;
    case 'B':
        break;
    default:
        Rf_error("Rsamtools internal: '%s' unhandled in _get_lkup", baseclass);
        break;
    }
    return lkup;
}

void _as_seqlevels(SEXP vec, SEXP lvls)
{
    
    _as_factor_SEXP(vec, lvls);
}

void _as_factor_SEXP(SEXP vec, SEXP lvls)
{
    SEXP cls = PROTECT(NEW_CHARACTER(1));
    SET_STRING_ELT(cls, 0, mkChar("factor"));
    SET_CLASS(vec, cls);
    SET_ATTR(vec, install("levels"), lvls);
    UNPROTECT(1);
}

KHASH_MAP_INIT_STR(str, int)

void _as_rname(SEXP vec, const char **lvls, const int n_lvls)
{
    /* map duplicate levels to first occurrence */
    SEXP levels = R_NilValue;
    int i, n = 0, ret, *map = NULL;
    khash_t(str) *level_map = kh_init(str);
    khiter_t iter;

    for (i = 0; i < n_lvls; ++i) {
        iter = kh_get(str, level_map, lvls[i]);
        if (iter == kh_end(level_map)) {
            iter = kh_put(str, level_map, lvls[i], &ret);
            kh_value(level_map, iter) = ++n;
        }
    }
    
    map = (int *) calloc(n_lvls, sizeof(int));
    levels = PROTECT(NEW_CHARACTER(n));
    for (i = 0; i < n_lvls; ++i) {
        iter = kh_get(str, level_map, lvls[i]);
        map[i] = kh_value(level_map, iter);
        SET_STRING_ELT(levels, map[i] - 1, mkChar(lvls[i]));
    }

    /* don't need to free keys (owned by calling function) */
    kh_destroy(str, level_map);

    int *v = INTEGER(vec);
    for (i = 0; i < Rf_length(vec); ++i)
        v[i] = v[i] == NA_INTEGER ? v[i] : map[v[i] - 1];

    _as_factor_SEXP(vec, levels);
    UNPROTECT(1);
}

void _as_factor(SEXP vec, const char **lvls, const int n_lvls)
{
    SEXP levels = PROTECT(NEW_CHARACTER(n_lvls));
    int i;
    for (i = 0; i < n_lvls; ++i)
        SET_STRING_ELT(levels, i, mkChar(lvls[i]));
    _as_factor_SEXP(vec, levels);
    UNPROTECT(1);
}

SEXP _as_XStringSet(const char **key, int len, const char *baseclass)
{
    char classname[40];         /* longest string should be "DNAStringSet" */

    if (snprintf(classname, sizeof(classname), "%sSet", baseclass)
        >= (int) sizeof(classname))
        error("Rsamtools internal error in _as_XStringSet(): "
              "'classname' buffer too small");

    SEXP lkup = _get_lkup(baseclass);
    int *lkup0, lkup_length = 0;
    if (R_NilValue == lkup) {
        lkup0 = NULL;
    } else {
        lkup0 = INTEGER(lkup);
        lkup_length = LENGTH(lkup);
    }

    SEXP width, ans;
    int i;
    PROTECT(width = NEW_INTEGER(len));
    for (i = 0; i < len; ++i)
        INTEGER(width)[i] = strlen(key[i]);
    PROTECT(ans = alloc_XRawList(classname, baseclass, width));

    XVectorList_holder holder;
    Chars_holder dest;
    holder = hold_XVectorList(ans);
    for (i = 0; i < len; ++i) {
        const char *seq = key[i];
        dest = get_elt_from_XRawList_holder(&holder, i);
        Ocopy_bytes_to_i1i2_with_lkup(0, dest.length - 1, (char *) dest.ptr,
                                      dest.length, seq, strlen(seq), lkup0,
                                      lkup_length);
    }

    UNPROTECT(2);
    return ans;
}

SEXP _as_PhredQuality(const char **key, int len)
{
    SEXP xstringset = PROTECT(_as_XStringSet(key, len, "BString"));

    SEXP s, t, nmspc, result;
    nmspc = PROTECT(_get_namespace("Rsamtools"));
    NEW_CALL(s, t, "PhredQuality", nmspc, 2);
    CSET_CDR(t, "x", xstringset);
    CEVAL_TO(s, nmspc, result);
    UNPROTECT(2);
    return result;
}

void _reverse(char *buf, int len)
{
    int i;
    char tmp;
    for (i = 0; i < floor(len / 2); ++i) {
        tmp = buf[len - i - 1];
        buf[len - i - 1] = buf[i];
        buf[i] = tmp;
    }
}

void _reverseComplement(char *buf, int len)
{
    static const int MAX_MAP = 256;
    static char map[256];
    static int init = 0;
    if (init == 0) {
        init = 1;
        for (int i = 0; i < MAX_MAP; ++i)
            map[i] = (char) i;
        map['A'] = 'T';
        map['C'] = 'G';
        map['G'] = 'C';
        map['T'] = 'A';
        map['a'] = 't';
        map['c'] = 'g';
        map['g'] = 'c';
        map['t'] = 'a';
        map['M'] = 'K';
        map['R'] = 'Y';
        map['Y'] = 'R';
        map['K'] = 'M';
        map['m'] = 'k';
        map['r'] = 'y';
        map['y'] = 'r';
        map['k'] = 'm';
        map['V'] = 'B';
        map['H'] = 'D';
        map['D'] = 'H';
        map['B'] = 'V';
        map['v'] = 'b';
        map['h'] = 'd';
        map['d'] = 'h';
        map['b'] = 'v';
    }
    _reverse(buf, len);
    for (int i = 0; i < len; ++i)
        buf[i] = map[(int) buf[i]];
}

char *_rtrim(char *s)
{
    int i = strlen(s) - 1;
    while (i >= 0) {
        if (s[i] != '\r')
            break;
        s[i--] = '\0';
    }
    return s;
}

/*
 * Copied from XVector/src/io_utils.c
 * Doesn't actually delete anything but returns the length that the 'buf' char
 * array would have had after deletion of the LF ("\n") or CRLF ("\r\n") chars
 * located at its end.
 * If 'buf_len' is -1, then 'buf' must be a C string (i.e. null-terminated).
 */
int _delete_trailing_LF_or_CRLF(const char *buf, int buf_len)
{
    if (buf_len == -1)
        buf_len = strlen(buf);
    if (buf_len == 0)
        return 0;
    if (buf[--buf_len] != '\n')
        return ++buf_len;
    if (buf_len == 0)
        return 0;
    if (buf[--buf_len] != '\r')
        return ++buf_len;
    return buf_len;
}

/* scan-related */

void _checkext(SEXP ext, SEXP tag, const char *lbl)
{
    if (EXTPTRSXP != TYPEOF(ext) || tag != R_ExternalPtrTag(ext))
        Rf_error("incorrect instance for '%s'", lbl);
}

void _checknames(SEXP filename, SEXP indexname, SEXP filemode)
{
    if (!IS_CHARACTER(filename) || LENGTH(filename) > 1)
        Rf_error("'filename' must be character(0) or character(1)");
    if (!IS_CHARACTER(indexname) || LENGTH(indexname) > 1)
        Rf_error("'indexname' must be character(0) or character(1)");
    if (!IS_CHARACTER(filemode) || LENGTH(filemode) != 1)
        Rf_error("'filemode' must be character(1)");
}

void _checkparams(SEXP regions, SEXP keepFlags, SEXP isSimpleCigar)
{
    const int MAX_CHRLEN = 1 << 29;	/* See samtools/bam_index.c */
    if (R_NilValue != regions) {
        if (!IS_VECTOR(regions) || LENGTH(regions) != 3)
            Rf_error("'regions' must be list(3) or NULL");
        if (!IS_CHARACTER(VECTOR_ELT(regions, 0)))
            Rf_error("internal: 'regions[1]' must be character()");
        if (!IS_INTEGER(VECTOR_ELT(regions, 1)))
            Rf_error("internal: 'regions[2]' must be integer()");
        if (!IS_INTEGER(VECTOR_ELT(regions, 2)))
            Rf_error("internal: 'regions[3]' must be integer()");
        if ((LENGTH(VECTOR_ELT(regions, 0)) != LENGTH(VECTOR_ELT(regions, 1))) ||
            (LENGTH(VECTOR_ELT(regions, 0)) != LENGTH(VECTOR_ELT(regions, 2))))
            Rf_error("internal: 'regions' elements must all be the same length");
        const int *end = INTEGER(VECTOR_ELT(regions, 2)),
            nrange = LENGTH(VECTOR_ELT(regions, 2));
        for (int irange = 0; irange < nrange; ++irange)
            if (end[irange] > MAX_CHRLEN)
                Rf_error("'end' must be <= %d", MAX_CHRLEN);
    }
    if (R_NilValue != keepFlags)
        if (!IS_INTEGER(keepFlags) || LENGTH(keepFlags) != 2)
            Rf_error("'keepFlags' must be integer(2) or NULL");
    if (R_NilValue != isSimpleCigar)
        if (!IS_LOGICAL(isSimpleCigar) || LENGTH(isSimpleCigar) != 1)
            Rf_error("'isSimpleCigar' must be logical(1) or NULL");
}

/* pairing */

static int check_x_or_y(SEXP qname, SEXP flag, SEXP rname, SEXP pos,
                        SEXP rnext, SEXP pnext, const char *what)
{
    int len;

    len = LENGTH(flag);
    if (qname != R_NilValue)
        if (!IS_CHARACTER(qname) || LENGTH(qname) != len)
            Rf_error("'%s_qname' must be NULL or a character vector "
                     "of the same length as '%s_flag'", what, what);
    if (!isFactor(rname) || LENGTH(rname) != len)
        Rf_error("'%s_rname' must be a factor "
                 "of the same length as '%s_flag'", what, what);
    if (!IS_INTEGER(pos) || LENGTH(pos) != len)
        Rf_error("'%s_pos' must be an integer vector "
                 "of the same length as '%s_flag'", what, what);
    if (!isFactor(rnext) || LENGTH(rnext) != len)
        Rf_error("'%s_rnext' must be a factor "
                 "of the same length as '%s_flag'", what, what);
    if (!IS_INTEGER(pnext) || LENGTH(pnext) != len)
        Rf_error("'%s_pnext' must be an integer vector "
                 "of the same length as '%s_flag'", what, what);
    return len;
}

/* 'x_qname' and 'y_qname' must be NULLs or 0-terminated strings. */
static int is_a_pair(const char *x_qname, int x_flag, int x_rname,
                     int x_pos, int x_rnext, int x_pnext,
                     const char *y_qname, int y_flag, int y_rname,
                     int y_pos, int y_rnext, int y_pnext)
{
    int x_is_not_first, x_is_not_last, x_is_proper, x_is_secondary,
        y_is_not_first, y_is_not_last, y_is_proper, y_is_secondary,
        tmp, ok;


    /* Flag bits 0x001, 0x004, and 0x008 must be 1, 0, and 0.
       Testing the 3 bits at once. */
    if (((x_flag & 0x00d) != 0x001) || ((y_flag & 0x00d) != 0x001))
        return 0;

    /* 'x' and 'y' must be the first or the last segment in the template. */
    x_is_not_first = (x_flag & 0x040) == 0;
    x_is_not_last = (x_flag & 0x080) == 0;
    y_is_not_first = (y_flag & 0x040) == 0;
    y_is_not_last = (y_flag & 0x080) == 0;
    if ((x_is_not_first == x_is_not_last) || (y_is_not_first == y_is_not_last))
        return 0;

    /* See top of the R/findMateAlignment.R file for the definition of tests
       (A), (B), (C), (D), (E), (F), and (G). */

    /* Test (A). */
    tmp = (x_qname != NULL) + (y_qname != NULL);
    if (tmp == 1)
        return 0;  /* 'x_qname' or 'y_qname' is NULL but not both */
    if (tmp == 2 && strcmp(x_qname, y_qname) != 0)
        return 0;  /* 'x_qname' and 'y_qname' differ */

    /* Test (B). */
    if ((x_rnext != y_rname) || (y_rnext != x_rname))
        return 0;

    /* Test (C). */
    if ((x_pnext != y_pos) || (y_pnext != x_pos))
        return 0;

    /* Test (D). */
    ok = (((x_flag & 0x020) == 0) == ((y_flag & 0x010) == 0)) &&
         (((y_flag & 0x020) == 0) == ((x_flag & 0x010) == 0));
    if (!ok)
        return 0;

    /* Test (E). */
    if (x_is_not_first == y_is_not_first)
        return 0;

    /* Test (F). */
    x_is_proper = (x_flag & 0x002) != 0;
    y_is_proper = (y_flag & 0x002) != 0;
    if (x_is_proper != y_is_proper)
        return 0;

    /* Test (G). */
    x_is_secondary = (x_flag & 0x100) != 0;
    y_is_secondary = (y_flag & 0x100) != 0;
    if (x_is_secondary != y_is_secondary)
        return 0;

    return 1;
}

/*
 * Parallel pairing of 2 vectors of alignments (i.e. BAM records).
 * The 2 input vectors 'x' and 'y' must have the same length.
 * 'x_rname', 'y_rname', 'x_rnext', and 'y_rnext' must all have exactly the
 * same levels in the same order. This is NOT checked.
 * Returns a logical vector of the same length as the input vectors.
 */
SEXP p_pairing(SEXP x_qname, SEXP x_flag, SEXP x_rname,
               SEXP x_pos, SEXP x_rnext, SEXP x_pnext,
               SEXP y_qname, SEXP y_flag, SEXP y_rname,
               SEXP y_pos, SEXP y_rnext, SEXP y_pnext)
{
    SEXP ans, x_qname_elt, y_qname_elt;
    int x_len, y_len, i,
        x_flag_elt, x_rname_elt, x_pos_elt, x_rnext_elt, x_pnext_elt,
        y_flag_elt, y_rname_elt, y_pos_elt, y_rnext_elt, y_pnext_elt;
    const char *x_qname_string, *y_qname_string;

    x_len = check_x_or_y(x_qname, x_flag, x_rname, x_pos,
                         x_rnext, x_pnext, "x");
    y_len = check_x_or_y(y_qname, y_flag, y_rname, y_pos,
                         y_rnext, y_pnext, "y");
    if (x_len != y_len)
        Rf_error("'x' and 'y' must have the same length");
    if ((x_qname == R_NilValue) != (y_qname == R_NilValue))
        Rf_error("both of 'x' and 'y' must either be NULL or not");
    if (x_qname == R_NilValue)
        x_qname_string = y_qname_string = NULL;

    PROTECT(ans = NEW_LOGICAL(x_len));
    for (i = 0; i < x_len; i++) {
        x_flag_elt = INTEGER(x_flag)[i];
        y_flag_elt = INTEGER(y_flag)[i];
        if (x_flag_elt == NA_INTEGER || y_flag_elt == NA_INTEGER) {
            UNPROTECT(1);
            Rf_error("'x_flag' or 'y_flag' contains NAs");
        }
        if (x_qname != R_NilValue) {
            x_qname_elt = STRING_ELT(x_qname, i);
            y_qname_elt = STRING_ELT(y_qname, i);
            if (x_qname_elt == NA_STRING || y_qname_elt == NA_STRING) {
                UNPROTECT(1);
                Rf_error("'x_qname' or 'y_qname' contains NAs");
            }
            x_qname_string = CHAR(x_qname_elt);
            y_qname_string = CHAR(y_qname_elt);
        }
        x_rname_elt = INTEGER(x_rname)[i];
        y_rname_elt = INTEGER(y_rname)[i];
        x_pos_elt = INTEGER(x_pos)[i];
        y_pos_elt = INTEGER(y_pos)[i];
        x_rnext_elt = INTEGER(x_rnext)[i];
        y_rnext_elt = INTEGER(y_rnext)[i];
        x_pnext_elt = INTEGER(x_pnext)[i];
        y_pnext_elt = INTEGER(y_pnext)[i];
        LOGICAL(ans)[i] = is_a_pair(x_qname_string, x_flag_elt, x_rname_elt,
                                    x_pos_elt, x_rnext_elt, x_pnext_elt,
                                    y_qname_string, y_flag_elt, y_rname_elt,
                                    y_pos_elt, y_rnext_elt, y_pnext_elt);
    }
    UNPROTECT(1);
    return ans;
}

/*
 * The input is a vector 'x' of alignments.
 * 'group_sizes 'must be a vector of positive integers that sum up to the
 * length of the input vectors. This is NOT checked.
 * 'x_rname' and 'x_rnext' must have exactly the same levels in the same order.
 * This is NOT checked.
 * Returns an integer vector of the same length as the input vectors.
 */
SEXP find_mate_within_groups(SEXP group_sizes,
                             SEXP x_flag, SEXP x_rname,
                             SEXP x_pos, SEXP x_rnext, SEXP x_pnext)
{
    SEXP ans;
    int x_len, *ans_p, ngroup, n, group_size, offset, i, i2, j, j2, ok,
        x_flag_elt, x_rname_elt, x_pos_elt, x_rnext_elt, x_pnext_elt,
        y_flag_elt, y_rname_elt, y_pos_elt, y_rnext_elt, y_pnext_elt,
        ans_i2_is_na, ans_j2_is_na;

    x_len = check_x_or_y(R_NilValue, x_flag, x_rname, x_pos,
                         x_rnext, x_pnext, "x");
    PROTECT(ans = NEW_INTEGER(x_len));
    for (i = 0, ans_p = INTEGER(ans); i < x_len; i++, ans_p++)
        *ans_p = NA_INTEGER;
    ngroup = LENGTH(group_sizes);
    for (n = offset = 0; n < ngroup; n++, offset += group_size) {
        group_size = INTEGER(group_sizes)[n];
        for (i = 1; i < group_size; i++) {
            i2 = i + offset;
            x_flag_elt = INTEGER(x_flag)[i2];
            if (x_flag_elt == NA_INTEGER) {
                UNPROTECT(1);
                Rf_error("'x_flag' contains NAs");
            }
            x_rname_elt = INTEGER(x_rname)[i2];
            x_pos_elt = INTEGER(x_pos)[i2];
            x_rnext_elt = INTEGER(x_rnext)[i2];
            x_pnext_elt = INTEGER(x_pnext)[i2];
            for (j = 0; j < i; j++) {
                j2 = j + offset;
                y_flag_elt = INTEGER(x_flag)[j2];
                if (y_flag_elt == NA_INTEGER) {
                    UNPROTECT(1);
                    Rf_error("'y_flag' contains NAs");
                }
                y_rname_elt = INTEGER(x_rname)[j2];
                y_pos_elt = INTEGER(x_pos)[j2];
                y_rnext_elt = INTEGER(x_rnext)[j2];
                y_pnext_elt = INTEGER(x_pnext)[j2];
                ok = is_a_pair(NULL, x_flag_elt, x_rname_elt,
                               x_pos_elt, x_rnext_elt, x_pnext_elt,
                               NULL, y_flag_elt, y_rname_elt,
                               y_pos_elt, y_rnext_elt, y_pnext_elt);
                if (!ok)
                    continue;
                ans_i2_is_na = INTEGER(ans)[i2] == NA_INTEGER;
                INTEGER(ans)[i2] = ans_i2_is_na ? j2 + 1 : 0;
                ans_j2_is_na = INTEGER(ans)[j2] == NA_INTEGER;
                INTEGER(ans)[j2] = ans_j2_is_na ? i2 + 1 : 0;
            }
        }
    }
    for (i = 0, ans_p = INTEGER(ans); i < x_len; i++, ans_p++) {
        if (*ans_p != NA_INTEGER && *ans_p != 0
         && INTEGER(ans)[*ans_p - 1] == 0)
            *ans_p = -(*ans_p);
    }
    UNPROTECT(1);
    return ans;
}
