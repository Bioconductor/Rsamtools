#include <stdlib.h>
#include <string.h>
#include "tabixfile.h"
#include "vcffile.h"

/* iterator to return null-terminated delimited fields */

struct it {
    char *str;
    char delim;
};

char *_it_next(struct it *it)
{
    char *curr = it->str;
    while ('\0' != *it->str && it->delim != *it->str)
        it->str++;
    if ('\0' != *it->str)
        *it->str++ = '\0';
    return curr;
}

char *_it_init(struct it *it, char *str, char delim)
{
    it->str = str;
    it->delim = delim;
    return _it_next(it);
}

/* parse a single 'vcf' line into vectors and matricies 'map' */

const struct fld_fmt {
    const char *name;
    SEXPTYPE type;
} FLD_FMT[] = {
    {
    "CHROM", STRSXP}, {
    "POS", INTSXP}, {
    "ID", STRSXP}, {
    "REF", STRSXP}, {
    "ALT", STRSXP}, {
    "QUAL", REALSXP}, {
    "FILTER", STRSXP}, {
    "INFO", STRSXP}, {
    "GENO", VECSXP}
};

const int N_FLDS = sizeof(FLD_FMT) / sizeof(FLD_FMT[0]);

SEXP _split_vcf(SEXP vcf, SEXP sample, SEXP imap, SEXP gmap)
{
    int i, j;
    const int
        vcf_n = Rf_length(vcf),
        samp_n = Rf_length(sample), 
        imap_n = Rf_length(imap),
        gmap_n = Rf_length(gmap),
        map_n = Rf_length(gmap);

    SEXP nms = GET_NAMES(gmap);
    SEXP inms = GET_NAMES(imap);
    int fmtidx, sampleidx;
    int *mapidx = (int *) R_alloc(sizeof(int), map_n),
        *imapidx = (int *) R_alloc(sizeof(int), imap_n);

    /* allocate result and first 7 fixed fields */
    SEXP result, info, geno, eltnms;
    PROTECT(result = Rf_allocVector(VECSXP, N_FLDS));
    for (i = 0; i < N_FLDS - 2; ++i)
        SET_VECTOR_ELT(result, i, Rf_allocVector(FLD_FMT[i].type, vcf_n));

    PROTECT(eltnms = Rf_allocVector(STRSXP, N_FLDS));
    for (i = 0; i < N_FLDS; ++i)
        SET_STRING_ELT(eltnms, i, mkChar(FLD_FMT[i].name));
    result = Rf_namesgets(result, eltnms);
    UNPROTECT(1);

    /* allocate info */
    PROTECT(info = Rf_allocVector(VECSXP, imap_n));
    for (j = 0; j < imap_n; ++j) {
        SEXPTYPE type = TYPEOF(VECTOR_ELT(imap, j));
        if (NILSXP == type) {
            SET_VECTOR_ELT(info, j, R_NilValue);
            continue;
        }
        SEXP elt = Rf_allocMatrix(type, vcf_n, 1);
        SET_VECTOR_ELT(info, j, elt);
        switch (type) {
        case INTSXP:
            for (i = 0; i < vcf_n; ++i)
                INTEGER(elt)[i] = R_NaInt;
            break;
        case REALSXP:
            for (i = 0; i < vcf_n; ++i)
                REAL(elt)[i] = R_NaReal;
            break;
        case STRSXP:
            for (i = 0; i < vcf_n; ++i)
                SET_STRING_ELT(elt, i, R_NaString);
            break;
        default:
            Rf_error("(internal) unhandled type '%s'", type2char(type));
        }
    }

    /* allocate geno, including matricies */
    PROTECT(geno = Rf_allocVector(VECSXP, gmap_n));
    PROTECT(eltnms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(eltnms, 0, R_NilValue);
    SET_VECTOR_ELT(eltnms, 1, sample);
    for (j = 0; j < map_n; ++j) {
        SEXPTYPE type = TYPEOF(VECTOR_ELT(gmap, j));
        if (NILSXP == type) {
            SET_VECTOR_ELT(geno, j, R_NilValue);
            continue;
        }
        SEXP elt = Rf_allocMatrix(type, vcf_n, samp_n);
        SET_VECTOR_ELT(geno, j, elt);
        switch (type) {
        case INTSXP:
            for (i = 0; i < vcf_n * samp_n; ++i)
                INTEGER(elt)[i] = R_NaInt;
            break;
        case REALSXP:
            for (i = 0; i < vcf_n * samp_n; ++i)
                REAL(elt)[i] = R_NaReal;
            break;
        case STRSXP:
            for (i = 0; i < vcf_n * samp_n; ++i)
                SET_STRING_ELT(elt, i, R_NaString);
            break;
        default:
            Rf_error("(internal) unhandled type '%s'", type2char(type));
        }
        elt = Rf_dimnamesgets(elt, eltnms);
    }
    UNPROTECT(1);

    /* parse each line */
    for (i = 0; i < vcf_n; i++) {
        struct it it0, it1, it2;
        char *record, *sample, *field, *ifld, *ikey, *ivlu, *fmt;

        record = strdup(CHAR(STRING_ELT(vcf, i)));

        /* first 7 'fixed' fields */
        for (field = _it_init(&it0, record, '\t'), j = 0;
             j < N_FLDS - 2; field = _it_next(&it0), ++j) {
            SEXP elt = VECTOR_ELT(result, j);
            switch (TYPEOF(elt)) {
            case INTSXP:
                INTEGER(elt)[i] = atoi(field);
                break;
            case REALSXP:
                REAL(elt)[i] = atof(field);
                break;
            case STRSXP:
                SET_STRING_ELT(elt, i, mkChar(field));
                break;
            default:
                Rf_error("(internal) unhandled fixed field type '%s'",
                         type2char(TYPEOF(elt)));
            }
        }

        /* 'INFO' field */
        ifld = field;
        for (field = _it_init(&it1, ifld, ';'), fmtidx = 0;
            '\0' != *field; field = _it_next(&it1), fmtidx++) {
            /* FIXME : allocate */
            char *txt = (char *) R_alloc(100, sizeof(char));
            strncpy(txt, field, 100);
            ikey = strtok(txt, "=");
            for (j = 0; j < imap_n; ++j) {
                if (0L == strcmp(ikey, CHAR(STRING_ELT(inms, j))))
                    break;
            }
            if (imap_n == j)
                Rf_error("record %d field %d INFO '%s' not found",
                         i + 1, fmtidx + 1, ikey);
            imapidx[fmtidx] = j;

            ivlu = strtok(NULL, "=");
            /* FIXME : set as "1" or "TRUE" */
            if (NULL == ivlu)
                ivlu = "1";
            SEXP matrix = VECTOR_ELT(info, imapidx[fmtidx]);
            int midx = i;
            switch (TYPEOF(matrix)) {
            case NILSXP:
                break;
            case INTSXP:
                INTEGER(matrix)[midx] = atoi(ivlu);
                break;
            case REALSXP:
                REAL(matrix)[midx] = atof(ivlu);
                break;
            case STRSXP:
                SET_STRING_ELT(matrix, midx, mkChar(ivlu));
                break;
            default:
                Rf_error("(internal) unhandled type '%s'",
                         type2char(TYPEOF(matrix)));
            }
        }
        field = _it_next(&it0);

        /* 'FORMAT' field */
        fmt = field;
        for (field = _it_init(&it2, fmt, ':'), fmtidx = 0;
             '\0' != *field; field = _it_next(&it2), fmtidx++) {
            for (j = 0; j < map_n; ++j) {
                if (0L == strcmp(field, CHAR(STRING_ELT(nms, j))))
                    break;
            }
            if (map_n == j)
                Rf_error("record %d field %d FORMAT '%s' not found",
                         i + 1, fmtidx + 1, field);
            mapidx[fmtidx] = j;
        }

        /* 'samples' field(s) */
        for (sample = _it_next(&it0), sampleidx = 0;
             '\0' != *sample; sample = _it_next(&it0), sampleidx++) {
            for (field = _it_init(&it2, sample, ':'), fmtidx = 0;
                 '\0' != *field; field = _it_next(&it2), fmtidx++) {
                SEXP matrix = VECTOR_ELT(geno, mapidx[fmtidx]);
                int midx = sampleidx * vcf_n + i;
                switch (TYPEOF(matrix)) {
                case NILSXP:
                    break;
                case INTSXP:
                    INTEGER(matrix)[midx] = atoi(field);
                    break;
                case REALSXP:
                    REAL(matrix)[midx] = atof(field);
                    break;
                case STRSXP:
                    SET_STRING_ELT(matrix, midx, mkChar(field));
                    break;
                default:
                    Rf_error("(internal) unhandled type '%s'",
                             type2char(TYPEOF(matrix)));
                }
            }
        }
        free(record);
    }

    /* remove NULL elements of info and geno  */
    for (i = 0, j = 0; i < Rf_length(geno); ++i) {
        if (R_NilValue != VECTOR_ELT(geno, i)) {
            SET_VECTOR_ELT(geno, j, VECTOR_ELT(geno, i));
            SET_STRING_ELT(nms, j, STRING_ELT(nms, i));
            j++; 
        }
    }
    PROTECT(nms = Rf_lengthgets(nms, j)); 
    PROTECT(geno = Rf_lengthgets(geno, j));
    geno = Rf_namesgets(geno, nms);
    UNPROTECT(2);

    for (i = 0, j = 0; i < Rf_length(info); ++i) {
        if (R_NilValue != VECTOR_ELT(info, i)) {
            SET_VECTOR_ELT(info, j, VECTOR_ELT(info, i));
            SET_STRING_ELT(inms, j, STRING_ELT(inms, i));
            j++; 
        }
    }
    PROTECT(inms = Rf_lengthgets(inms, j)); 
    PROTECT(info = Rf_lengthgets(info, j));
    info = Rf_namesgets(info, inms);
    UNPROTECT(2);

    SET_VECTOR_ELT(result, N_FLDS - 2, info);
    SET_VECTOR_ELT(result, N_FLDS - 1, geno);

    UNPROTECT(3);
    return result;
}

SEXP scan_vcf(SEXP ext, SEXP space, SEXP yieldSize, SEXP sample, 
              SEXP imap, SEXP gmap)
{
    SEXP tbx = PROTECT(scan_tabix(ext, space, yieldSize));
 
    for (int i = 0; i < Rf_length(tbx); ++i) {
        SEXP result = _split_vcf(VECTOR_ELT(tbx, i), sample, imap, gmap);
        SET_VECTOR_ELT(tbx, i, result);
    }

    UNPROTECT(1);
    return tbx;
}

SEXP scan_vcf_connection(SEXP txt, SEXP sample, SEXP imap, SEXP gmap)
{
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 1));

    SET_VECTOR_ELT(result, 0, _split_vcf(txt, sample, imap, gmap));

    UNPROTECT(1);
    return result;
}
