noinst_LTLIBRARIES = libsphinxutil.la

libsphinxutil_la_SOURCES = \
        bitvec.c \
	bio.c \
        case.c \
        ckd_alloc.c \
	dtoa.c \
	listelem_alloc.c \
        cmd_ln.c \
        err.c \
        filename.c \
	garray.c \
	genrand.c \
        glist.c \
	gq.c \
        hash_table.c \
        heap.c \
	huff_code.c \
        info.c \
	logmath.c \
	mmio.c \
        pio.c \
        profile.c \
	sbthread.c \
        strfuncs.c \
	sync_array.c \
	utf8.c

noinst_HEADERS = garray_internal.h

INCLUDES = -I$(top_srcdir)/include \
	 -I$(top_srcdir)/include/sphinxbase \
	-I$(top_builddir)/include 

LIBOBJS = @LIBOBJS@
