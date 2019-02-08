#include "tagfilter.h"

C_TAGFILTER _tagFilter_as_C_types(SEXP tl) {
    if(LENGTH(tl) == 0)
        return NULL;
    C_TAGFILTER tagfilt = Calloc(1, _C_TAGFILTER);
    SEXP nmssxp = Rf_getAttrib(tl, R_NamesSymbol);
    /* convert tagnames */
    int len = LENGTH(nmssxp);
    tagfilt->len = len;
    tagfilt->tagnames = (const char**) Calloc(len, char*);
    for(int i = 0; i < len; ++i) {
        tagfilt->tagnames[i] = CHAR(STRING_ELT(nmssxp, i));
    }

    /* allocate for elements */
    tagfilt->elts = Calloc(len, _TAGFILTER_ELT);

    for(int i = 0; i < len; ++i) {
        SEXP sxp_elt = VECTOR_ELT(tl, i);
        int elt_len = LENGTH(sxp_elt);
        if(elt_len < 1)
            Rf_error("elements of tag filter list must have non-zero length");
        switch(TYPEOF(sxp_elt)) {
        case INTSXP: /* INTEGER */
            tagfilt->elts[i].len = elt_len;
            tagfilt->elts[i].type = TAGFILT_T_INT;
            tagfilt->elts[i].ptr = (int*) INTEGER(sxp_elt);
            break;
        case STRSXP: /* CHAR */
            tagfilt->elts[i].len = elt_len;
            tagfilt->elts[i].type = TAGFILT_T_STRING;
            tagfilt->elts[i].ptr = (const char**) Calloc(elt_len, char*);
            for(int j = 0; j < elt_len; ++j) {
                ((const char**) tagfilt->elts[i].ptr)[j] =
                    CHAR(STRING_ELT(sxp_elt, j));
            }
            break;
        default:
            Rf_error("unpermitted tag filter input type '%s'",
                     Rf_type2char((SEXPTYPE) TYPEOF(sxp_elt)));
        }
    }
    return tagfilt;
}

void _Free_C_TAGFILTER(C_TAGFILTER ctf) {
    if(ctf) {
        Free(ctf->tagnames);
        if(ctf->elts) {
            for(int i = 0; i < ctf->len; ++i) {
                if(ctf->elts[i].type == TAGFILT_T_STRING)
                    Free(ctf->elts[i].ptr);
            }
            Free(ctf->elts);
        }
        Free(ctf);
    }
}

static const char* const TagFilterType_str[]  = { "INTERNAL_ERROR: UNSET",
                                           "integer()", "character()" };

static const char auxtype[] = "cCsSiIfdAZHB";
static const char inttype[] = "cCsSiI";
static const char* const auxtype_str[] = {
    "integer", "integer", "integer", "integer", "integer", "integer",
    "float", "double", "single printable character", "string",
    "hex array", "array" };

static void _typemismatch_error(const char *tagname, uint8_t *aux,
                                TagFilterType tft, const char* val_as_string,
                                int irec) {
    static const char* msg =
        "tag '%s' type ('%s') does not match tagFilter type\n"
        "    BAM read tag:   %s:%c:%s\n"
        "    tagFilter type: %s\n"
        "    Record number:  %d";
    int i = strchr(auxtype, aux[0]) - auxtype;
    const char *printable_typename = auxtype_str[i];
    const char tagtypechar = strchr(inttype, aux[0]) ? 'i' : aux[0];
    Rf_error(msg, tagname, printable_typename, tagname, tagtypechar,
             val_as_string, TagFilterType_str[tft], irec);
}

static void _typeunsupported_error(const char *tagname, uint8_t *aux,
    const char* val_as_string, int irec) {
    static const char* msg =
        "tag '%s' type ('%s') unsupported by tagFilter\n"
        "    BAM read tag:  %s:%c:%s\n"
        "    Record number: %d";
    int i = strchr(auxtype, aux[0]) - auxtype;
    const char *printable_typename = auxtype_str[i];
    const char tagtypechar = strchr(inttype, aux[0]) ? 'i' : aux[0];
    Rf_error(msg, tagname, printable_typename, tagname, tagtypechar, val_as_string, irec);
}

int _tagfilter(const bam1_t * bam, C_TAGFILTER tagfilter, int irec)
{
    int len = tagfilter->len;
    for(int i = 0; i < len; ++i) {
        const char* tagname = tagfilter->tagnames[i];
        void *list_elt = tagfilter->elts[i].ptr;
        uint8_t *aux = bam_aux_get(bam, tagname);
        if(aux == NULL)
            return 0;

        int idx;
        /* one bamf* for each case in switch */
        int bamfi;
        float bamff;
        double bamfd;
        char bamfA;
        char* bamfZ;
        char val_as_string[51];
        switch(aux[0]) {
        case 'c':
        case 'C':
        case 's':
        case 'S':
        case 'i':
        case 'I': /* INTEGER */
            bamfi = bam_aux2i(aux);
            if(tagfilter->elts[i].type != TAGFILT_T_INT) {
                snprintf(val_as_string, 51, "%d", bamfi);
                _typemismatch_error(tagname, aux, tagfilter->elts[i].type,
                                    val_as_string, irec);
            }
            for(idx = 0; idx < tagfilter->elts[i].len; ++idx) {
                if(bamfi == ((int*) list_elt)[idx])
                    break;
            }
            if(idx == tagfilter->elts[i].len)
                return 0;
            break;
        case 'f': /* REAL */
            bamff = (float) bam_aux2f(aux);
            snprintf(val_as_string, 51, "%f", bamff);
            _typeunsupported_error(tagname, aux, val_as_string, irec);
            break;
        case 'd': /* REAL */
            bamfd = bam_aux2f(aux);
            snprintf(val_as_string, 51, "%f", bamfd);
            _typeunsupported_error(tagname, aux, val_as_string, irec);
            break;
        case 'A': /* STRSXP */
            bamfA = bam_aux2A(aux);
            if(tagfilter->elts[i].type != TAGFILT_T_STRING ||
               strlen(*((const char**) list_elt)) != 1) {
                snprintf(val_as_string, 51, "%c", bamfA);
                _typemismatch_error(tagname, aux, tagfilter->elts[i].type,
                                    val_as_string, irec);
            }
            for(idx = 0; idx < tagfilter->elts[i].len; ++idx) {
                if(bamfA == *((const char**) list_elt)[idx])
                    break;
            }
            if(idx == tagfilter->elts[i].len)
                return 0;
            break;
        case 'Z': /* STRSXP */
            bamfZ = bam_aux2Z(aux);
            if(tagfilter->elts[i].type != TAGFILT_T_STRING) {
                snprintf(val_as_string, 51, "%s", bamfZ);
                _typemismatch_error(tagname, aux, tagfilter->elts[i].type,
                                    val_as_string, irec);
            }
            for(idx = 0; idx < tagfilter->elts[i].len; ++idx) {
                if(!strcmp(bamfZ, ((const char**) list_elt)[idx]))
                    break;
            }
            if(idx == tagfilter->elts[i].len)
                return 0;
            break;
        case 'H': /* RAWSXP */
            /* Can actually use bam_aux2Z(aux) to get a
             * null-terminated char array */
            bamfZ = bam_aux2Z(aux);
            snprintf(val_as_string, 51, "%s", bamfZ);
            _typeunsupported_error(tagname, aux, val_as_string, irec);
            break;
        case 'B':
            _typeunsupported_error(tagname, aux, "[unknown]", irec);
            break;
        default:
            Rf_error("unknown tag type '%c', record %d", aux[0], irec);
            break;
        }
    }

    return 1;
}
