# over-written except for Windows
SAMTOOLS_PATH=\
    $(shell echo 'cat(system.file("usrlib", .Platform[["r_arch"]],\
                         package="Rsamtools", mustWork=TRUE))' |\
                 "${R_HOME}/bin/R" --vanilla --slave)
PKG_LIBS+=-L"$(SAMTOOLS_PATH)" -lbam -lbcf $(ZLIB_LIBS) -lws2_32
PKG_CPPFLAGS+=-D_USE_KNETFILE -D_FILE_OFFSET_BITS=64 \
	-D_LARGEFILE64_SOURCE
