check_PROGRAMS = \
	test_fsg_read \
	test_fsg_jsgf \
	test_fsg_write_fsm

TESTS = $(check_PROGRAMS)

INCLUDES = -I$(top_srcdir) \
	-DLMDIR=\"${srcdir}\"

LDADD = ${top_builddir}/libmultisphinx.la

noinst_HEADERS = test_macros.h

EXTRA_DIST = goforward.fsg polite.gram
