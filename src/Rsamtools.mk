# over-written except for Windows
SAMTOOLS_PATH=\
    `echo 'cat(system.file("usrlib", .Platform[["r_arch"]],\
                         package="Rsamtools", mustWork=TRUE))' |\
                 "${R_HOME}/bin/R" --vanilla --slave`
PKG_LIBS+=-L"$(SAMTOOLS_PATH)" -lbam -lbcf -ltabix -lws2_32
PKG_CPPFLAGS+=-D_USE_KNETFILE -D_FILE_OFFSET_BITS=64 \
	-D_LARGEFILE64_SOURCE

# zlib
ZLIB_CFLAGS+=`echo 'zlibbioc::pkgconfig("PKG_CFLAGS")'|\
    "${R_HOME}/bin/R" --vanilla --slave`
PKG_LIBS+=`echo 'zlibbioc::pkgconfig("PKG_LIBS_shared")' |\
    "${R_HOME}/bin/R" --vanilla --slave`
%.o: %.c
	$(CC) $(ZLIB_CFLAGS) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@

