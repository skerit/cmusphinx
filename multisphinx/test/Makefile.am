check_PROGRAMS =				\
	test_arc_buffer				\
	test_bptbl				\
	test_build_lattice			\
	test_forward_backward			\
	test_htk_lattice			\
	test_latgen				\
	test_ngram_trie				\
	test_ngram_trie_iteration		\
	test_nodeid_map				\
	test_partial_backward			\
	test_search_factory			\
	test_state_align		\
	test_partial_results

TESTS = $(check_PROGRAMS)

LDADD = $(top_builddir)/libmultisphinx.la
INCLUDES = -I$(top_srcdir) -DTESTDATADIR=\"$(abs_srcdir)/data\"
