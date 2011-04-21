check_PROGRAMS = \
	test_bptbl \
	test_fwdtree \
	test_fwdtree_fwdflat \
	test_fwdtree_fwdflat_homos \
	test_latgen

TESTS = $(check_PROGRAMS)

INCLUDES = \
	-I$(top_srcdir) \
	-DTESTDATADIR=\"$(top_srcdir)/test/data\" \
	-DMODELDIR=\"$(top_srcdir)/model\"

LDADD = ${top_builddir}/libmultisphinx.la

noinst_HEADERS = test_macros.h

CLEANFILES = *.log
