#include "test_macros.h"

#include <multisphinx/ngram_trie.h>
#include <multisphinx/pocketsphinx_internal.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int
main(int argc, char *argv[])
{
	ngram_trie_node_t *root, *ug1, *ug2, *bg1, *bg2, *bg3, *tg;
	ngram_trie_iter_t *itor;
	ngram_trie_t *t;
	cmd_ln_t *config;
	logmath_t *lmath;
	bin_mdef_t *mdef;
	dict_t *dict;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-dict", TESTDATADIR "/bn10000.homos.dic",
			     NULL);
	ps_init_defaults(config);

	lmath = logmath_init(cmd_ln_float32_r(config, "-logbase"),
			     0, FALSE);
	mdef = bin_mdef_read(config, cmd_ln_str_r(config, "-mdef"));
	dict = dict_init(config, mdef);
	t = ngram_trie_init(dict, lmath);

	root = ngram_trie_root(t);
	ug1 = ngram_trie_add_successor(t, root, 1);
	ug2 = ngram_trie_add_successor(t, root, 2);
	bg1 = ngram_trie_add_successor(t, ug1, 3);
	bg2 = ngram_trie_add_successor(t, ug1, 4);
	bg3 = ngram_trie_add_successor(t, ug2, 5);
	tg = ngram_trie_add_successor(t, bg3, 6);

	TEST_ASSERT(3 == ngram_trie_n(t));

	itor = ngram_trie_ngrams(t, 1);
	TEST_ASSERT(ngram_trie_iter_get(itor) == ug1);
	itor = ngram_trie_iter_next(itor);
	TEST_ASSERT(ngram_trie_iter_get(itor) == ug2);
	itor = ngram_trie_iter_next(itor);
	TEST_ASSERT(itor == NULL);

	itor = ngram_trie_ngrams(t, 2);
	TEST_ASSERT(ngram_trie_iter_get(itor) == bg1);
	itor = ngram_trie_iter_next(itor);
	TEST_ASSERT(ngram_trie_iter_get(itor) == bg2);
	itor = ngram_trie_iter_next(itor);
	TEST_ASSERT(ngram_trie_iter_get(itor) == bg3);
	itor = ngram_trie_iter_next(itor);
	TEST_ASSERT(itor == NULL);

	itor = ngram_trie_ngrams(t, 3);
	TEST_ASSERT(ngram_trie_iter_get(itor) == tg);
	itor = ngram_trie_iter_next(itor);
	TEST_ASSERT(itor == NULL);

	ngram_trie_free(t);
	return 0;
}
