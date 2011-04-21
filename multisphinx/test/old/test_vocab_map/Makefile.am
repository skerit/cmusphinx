check_PROGRAMS = \
	test_vocab_map

TESTS = $(check_PROGRAMS)

INCLUDES = -I$(top_srcdir) \
	-DTESTDATADIR=\"$(top_srcdir)/test/data\"

LDADD = ${top_builddir}/libmultisphinx.la

noinst_HEADERS = test_macros.h

CLEANFILES = *.log
