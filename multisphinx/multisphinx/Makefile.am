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
	ps_alignment.c				\
	ps_search.c				\
	ptm_mgau.c				\
	sendump.c				\
	s2_semi_mgau.c				\
	tmat.c					\
	vector.c				\
	vocab_map.c				\
	pocketsphinx.c

noinst_HEADERS =				\
	pocketsphinx_internal.h			\
	hmm.h					\
	mdef.h					\
	ms_gauden.h				\
	ms_lattice.h				\
	ms_mgau.h				\
	ms_senone.h				\
	nodeid_map.h				\
	posixwin32.h				\
	ptm_mgau.h				\
	sendump.h				\
	s2_semi_mgau.h				\
	tied_mgau_common.h			\
	tmat.h					\
	vector.h

# Need this to avoid pulling in system includes by accident
INCLUDES = -I$(top_srcdir)
