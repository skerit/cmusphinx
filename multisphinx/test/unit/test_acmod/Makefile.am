check_PROGRAMS = \
	test_sync_array \
	test_featbuf \
	test_acmod

TESTS = $(check_PROGRAMS)

INCLUDES = \
	-I$(top_srcdir) \
	-DTESTDATADIR=\"$(top_srcdir)/test/data\" \
	-DMODELDIR=\"$(top_srcdir)/model\"

LDADD = ${top_builddir}/libmultisphinx.la

noinst_HEADERS = test_macros.h

CLEANFILES = *.log
