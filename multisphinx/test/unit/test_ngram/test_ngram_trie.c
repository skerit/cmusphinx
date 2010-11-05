#include "ngram_trie.h"

#include "test_macros.h"
#include "pocketsphinx_internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int
main(int argc, char *argv[])
{
	ngram_trie_t *t;
	dict_t *dict;
	bin_mdef_t *mdef;
	logmath_t *lmath;
	cmd_ln_t *config;

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

	ngram_trie_free(t);
	dict_free(dict);
	logmath_free(lmath);
	bin_mdef_free(mdef);

	return 0;
}
