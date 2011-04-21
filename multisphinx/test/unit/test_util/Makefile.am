check_PROGRAMS = \
	test_fopen \
	test_huff_code \
	test_bit_encode \
	test_build_directory \
	test_heap \
	test_garray \
	test_gq \
	test_filename

TESTS = $(check_PROGRAMS)

INCLUDES = -I$(top_srcdir) \
	-DLMDIR=\"${top_srcdir}/test/unit/test_ngram\"

LDADD = ${top_builddir}/libmultisphinx.la

noinst_HEADERS = test_macros.h

CLEANFILES = *.log *.out
