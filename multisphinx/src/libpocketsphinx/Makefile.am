ARCH   = `uname -m | sed s/i.86/i386/`_`uname -s | tr A-Z a-z`

lib_LTLIBRARIES = libpocketsphinx.la

# Confusing libtool version info means:
#
# This is API version 3, software release 0, which only supports API 3
libpocketsphinx_la_LDFLAGS = -version-info 3:0:0

libpocketsphinx_la_LIBADD = $(top_builddir)/src/libsphinxbase/libsphinxbase.la

libpocketsphinx_la_SOURCES =			\
	acmod.c					\
	bin_mdef.c				\
	bptbl.c					\
	dict.c					\
	dict2pid.c				\
	featbuf.c				\
	fwdtree_search.c			\
	fwdflat_search.c			\
	hmm.c					\
	mdef.c					\
	ms_gauden.c				\
	ms_mgau.c				\
	ms_senone.c				\
	ps_alignment.c				\
	ps_search.c				\
	ptm_mgau.c				\
	sendump.c				\
	s2_semi_mgau.c				\
	state_align_search.c			\
	tmat.c					\
	vector.c				\
	pocketsphinx.c

noinst_HEADERS =				\
	pocketsphinx_internal.h			\
	acmod.h					\
	bptbl.h					\
	bin_mdef.h				\
	dict.h					\
	dict2pid.h				\
	featbuf.h				\
	fwdtree_search.h			\
	fwdflat_search.h			\
	hmm.h					\
	mdef.h					\
	ms_gauden.h				\
	ms_mgau.h				\
	ms_senone.h				\
	ps_alignment.h				\
	ps_search.h				\
	posixwin32.h				\
	ptm_mgau.h				\
	sendump.h				\
	s2_semi_mgau.h				\
	s3types.h				\
	tied_mgau_common.h			\
	tmat.h					\
	vector.h

INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/include \
	   -DMODELDIR=\"${prefix}/share/pocketsphinx/model\"

