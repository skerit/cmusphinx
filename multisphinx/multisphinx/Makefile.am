noinst_LTLIBRARIES = libmultisphinx.la
libmultisphinx_la_LDFLAGS = -version-info 3:0:0
libmultisphinx_la_SOURCES =			\
	acmod.c					\
	arc_buffer.c				\
	bin_mdef.c				\
	bptbl.c					\
	dict.c					\
	dict2pid.c				\
	featbuf.c				\
	fwdtree_search.c			\
	fwdflat_search.c			\
	hmm.c					\
	latgen_search.c				\
	mdef.c					\
	ms_gauden.c				\
	ms_lattice.c				\
	ms_mgau.c				\
	ms_senone.c				\
	nodeid_map.c				\
	ngram_trie.c				\
	alignment.c				\
	search.c				\
	ptm_mgau.c				\
	sendump.c				\
	s2_semi_mgau.c				\
	state_align_search.c	\
	tmat.c					\
	vector.c				\
	vocab_map.c				\
	search_factory.c

pkginclude_HEADERS = \
	acmod.h					\
	arc_buffer.h				\
	bin_mdef.h				\
	bptbl.h					\
	cmdln_macro.h				\
	dict.h					\
	dict2pid.h				\
	featbuf.h				\
	fwdflat_search.h			\
	fwdtree_search.h			\
	hmm.h					\
	latgen_search.h				\
	mdef.h					\
	ms_gauden.h				\
	ms_lattice.h				\
	ms_mgau.h				\
	ms_senone.h				\
	ngram_trie.h				\
	nodeid_map.h				\
	win32_export.h			\
	posixwin32.h				\
	alignment.h				\
	search.h				\
	search_internal.h			\
	ptm_mgau.h				\
	s2_semi_mgau.h				\
	s3types.h				\
	sendump.h				\
	state_align_search.h			\
	tied_mgau_common.h			\
	tmat.h					\
	vector.h				\
	vocab_map.h				\
	search_factory.h

# Need this to avoid pulling in system includes by accident
INCLUDES = -I$(top_srcdir)
