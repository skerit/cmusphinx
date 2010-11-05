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
	FILE *arpafh;
	int32 prob;
	int n_used;

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
	arpafh = fopen(TESTDATADIR "/bn10000.3g.arpa", "r");
	ngram_trie_read_arpa(t, arpafh);
	fclose(arpafh);

	/* Test 1, 2, 3-gram probs without backoff. */
	prob = ngram_trie_prob(t, &n_used, "THREE", "POINT", "ZERO", NULL);
	printf("P(ZERO POINT THREE) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -25776);
	prob = ngram_trie_prob(t, &n_used, "THREE", "POINT", NULL);
	printf("P(POINT THREE) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -38960);
	prob = ngram_trie_prob(t, &n_used, "THREE", NULL);
	printf("P(THREE) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -69328);

	/* Test 3-gram probs with backoff. */
	/* Backoff to 2-gram POINT FOUR + alpha(ZERO POINT) */
	prob = ngram_trie_prob(t, &n_used, "FOUR", "POINT", "ZERO", NULL);
	printf("P(ZERO POINT FOUR) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -35600);
	/* Backoff to 2-gram SIX FOUR + alpha(ZERO) */
	prob = ngram_trie_prob(t, &n_used, "FOUR", "SIX", "ZERO", NULL);
	printf("P(ZERO SIX FOUR) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -56608);
	/* Backoff to 1-gram FOUR + alpha(ZERO SEVEN) */
	prob = ngram_trie_prob(t, &n_used, "FOUR", "SEVEN", "ZERO", NULL);
	printf("P(ZERO SEVEN FOUR) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -76496);

	arpafh = fopen("tmp.bn10000.3g.arpa", "w");
	ngram_trie_write_arpa(t, arpafh);
	fclose(arpafh);
	ngram_trie_free(t);

	t = ngram_trie_init(dict, lmath);
	arpafh = fopen("tmp.bn10000.3g.arpa", "r");
	ngram_trie_read_arpa(t, arpafh);
	fclose(arpafh);

	/* Test 1, 2, 3-gram probs without backoff. */
	prob = ngram_trie_prob(t, &n_used, "THREE", "POINT", "ZERO", NULL);
	printf("P(ZERO POINT THREE) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -25776);
	prob = ngram_trie_prob(t, &n_used, "THREE", "POINT", NULL);
	printf("P(POINT THREE) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -38960);
	prob = ngram_trie_prob(t, &n_used, "THREE", NULL);
	printf("P(THREE) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -69328);

	/* Test 3-gram probs with backoff. */
	/* Backoff to 2-gram POINT FOUR + alpha(ZERO POINT) */
	prob = ngram_trie_prob(t, &n_used, "FOUR", "POINT", "ZERO", NULL);
	printf("P(ZERO POINT FOUR) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -35600);
	/* Backoff to 2-gram SIX FOUR + alpha(ZERO) */
	prob = ngram_trie_prob(t, &n_used, "FOUR", "SIX", "ZERO", NULL);
	printf("P(ZERO SIX FOUR) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -56608);
	/* Backoff to 1-gram FOUR + alpha(ZERO SEVEN) */
	prob = ngram_trie_prob(t, &n_used, "FOUR", "SEVEN", "ZERO", NULL);
	printf("P(ZERO SEVEN FOUR) = %d = %g = %f\n",
	       prob, logmath_exp(lmath, prob), logmath_log_to_log10(lmath, prob));
	TEST_EQUAL_LOG(prob, -76496);

	ngram_trie_free(t);
	dict_free(dict);
	logmath_free(lmath);
	bin_mdef_free(mdef);

	return 0;
}
