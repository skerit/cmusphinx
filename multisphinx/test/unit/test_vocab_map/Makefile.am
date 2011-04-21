check_PROGRAMS = \
	test_vocab_map

TESTS = $(check_PROGRAMS)

INCLUDES = \
	-I$(top_srcdir)/include/sphinxbase \
	-I$(top_srcdir)/include/multisphinx \
	-I$(top_srcdir)/src/libpocketsphinx \
	-I$(top_builddir)/include \
	-DTESTDATADIR=\"$(top_srcdir)/test/data\"

LDADD = ${top_builddir}/src/libsphinxbase/libsphinxbase.la \
	 ${top_builddir}/src/libpocketsphinx/libpocketsphinx.la

noinst_HEADERS = test_macros.h

CLEANFILES = *.log
