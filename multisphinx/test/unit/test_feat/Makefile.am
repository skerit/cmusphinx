check_PROGRAMS = test_feat test_feat_live test_feat_fe test_subvq
noinst_HEADERS = test_macros.h

INCLUDES = -I$(top_srcdir) \
	-DTESTDATADIR=\"$(top_srcdir)/test/data\"


LDADD = ${top_builddir}/libmultisphinx.la

TESTS = _test_feat.test test_feat_live test_feat_fe test_subvq
EXTRA_DIST = _test_feat.res _test_feat.test
CLEANFILES = *.out
