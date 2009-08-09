AC_INIT("DESCRIPTION")

AC_CHECK_LIB([z], [gzeof],,[AC_MSG_ERROR([need libz with gzeof])])

dnl for SAMTOOLS
AC_CHECK_FUNCS([bzero memset pow sqrt strcasecmp strdup strstr strtol])
AC_CHECK_HEADERS([zlib.h fcntl.h malloc.h stddef.h])
AC_CHECK_TYPES([ptrdiff_t])
AC_C_CONST
AC_C_INLINE
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_HEADER_STDBOOL
AC_TYPE_OFF_T
AC_TYPE_SIZE_T

AC_OUTPUT(src/Makevars)
