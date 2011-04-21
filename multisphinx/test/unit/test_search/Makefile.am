check_PROGRAMS = \
	test_bptbl \
	test_fwdtree \
	test_fwdtree_fwdflat \
	test_fwdtree_fwdflat_homos \
	test_latgen

TESTS = $(check_PROGRAMS)

INCLUDES = \
	-I$(top_srcdir)/include/sphinxbase \
	-I$(top_srcdir)/include/multisphinx \
	-I$(top_srcdir)/src/libpocketsphinx \
	-I$(top_builddir)/include \
	-DTESTDATADIR=\"$(top_srcdir)/test/data\" \
	-DMODELDIR=\"$(top_srcdir)/model\"

LDADD = ${top_builddir}/src/libsphinxbase/libsphinxbase.la \
	 ${top_builddir}/src/libpocketsphinx/libpocketsphinx.la

noinst_HEADERS = test_macros.h

CLEANFILES = *.log
