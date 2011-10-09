#include "encode.h"
#include "Biostrings_interface.h"

unsigned char _bDecode(char);
unsigned char _dnaDecode(char);
unsigned char _rnaDecode(char);

char _bEncode(char c)
{
    return c;
}

#define _dnaEncode DNAencode;
#define _rnaEncode RNAencode;

ENCODE_FUNC _encoder(const char *base)
{
    ENCODE_FUNC encode = NULL;
    if (strcmp(base, "DNAString") == 0) {
        encode = _dnaEncode;
    } else if (strcmp(base, "RNAString") == 0) {
        encode = _rnaEncode;
    } else if (strcmp(base, "BString") == 0) {
        encode = _bEncode;
    } else {
        Rf_error("internal: unknown '_encoder' class '%s'", base);
    }
    return encode;
}
