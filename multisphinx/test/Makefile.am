check_PROGRAMS =				\
	test_arc_buffer				\
	test_bptbl				\
	test_build_lattice			\
	test_forward_backward			\
	test_fwdflat				\
	test_fwdtree				\
	test_fwdtree_fwdflat			\
	test_fwdtree_fwdflat_homos		\
	test_htk_lattice			\
	test_latgen				\
	test_ngram_trie				\
	test_ngram_trie_iteration		\
	test_nodeid_map				\
	test_partial_backward

TESTS = $(check_PROGRAMS)

LDADD = $(top_builddir)/libmultisphinx.la
INCLUDES = -I$(top_srcdir) -DTESTDATADIR=\"$(srcdir)/data\"

