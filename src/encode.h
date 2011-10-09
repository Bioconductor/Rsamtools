#ifndef _RENCODE_H_
#define _RENCODE_H_

typedef char (*ENCODE_FUNC) (char);	/* DNAencode, RNAencode */
ENCODE_FUNC _encoder(const char *baseclass);	/* RNAString, DNAString, BString */

#endif
