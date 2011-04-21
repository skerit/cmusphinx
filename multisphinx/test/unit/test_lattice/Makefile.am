check_PROGRAMS = \
	test_nodeid_map \
	test_htk_lattice \
	test_build_lattice \
	test_forward_backward \
	test_partial_backward

TESTS = $(check_PROGRAMS)

INCLUDES = -I$(top_srcdir) \
	-DTESTDATADIR=\"$(top_srcdir)/test/data\" \
	-DMODELDIR=\"$(top_srcdir)/model\"

LDADD = ${top_builddir}/libmultisphinx.la

noinst_HEADERS = test_macros.h

CLEANFILES = *.log test_htk_lattice.slf test_htk_lattice.dot
