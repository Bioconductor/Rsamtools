SAMTOOLS_PATH=\
    `echo 'cat(system.file("usrlib", .Platform[["r_arch"]],\
                         package="Rsamtools", mustWork=TRUE))' |\
                 "${R_HOME}/bin/R" --vanilla --slave`
SAMTOOLS_LIBS="$(SAMTOOLS_PATH)/libbam.a" "$(SAMTOOLS_PATH)/libbcf.a"\
	"$(SAMTOOLS_PATH)/libtabix.a" -lz -pthread
SAMTOOLS_CPPFLAGS=-D_USE_KNETFILE -DBGZF_CACHE -D_FILE_OFFSET_BITS=64 \
	-D_LARGEFILE64_SOURCE
