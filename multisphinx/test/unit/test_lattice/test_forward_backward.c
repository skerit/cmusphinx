#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <multisphinx/ms_lattice.h>
#include <sphinxbase/ngram_model.h>

#include "test_macros.h"

int
main(int argc, char *argv[])
{
	ms_lattice_t *l;
	ngram_model_t *lm;
	logmath_t *lmath;
	FILE *fh;

	lmath = logmath_init(1.0001, 0, FALSE);
	TEST_ASSERT(lm = ngram_model_read(NULL, TESTDATADIR "/bn10000.3g.arpa",
					  NGRAM_ARPA, lmath));
	TEST_ASSERT(l = ms_lattice_init(lmath, NULL));

	/* Read in the test lattice. */
	TEST_ASSERT(fh = fopen(TESTDATADIR "/050c0103.slf", "r"));
	TEST_ASSERT(0 == ms_lattice_read_htk(l, fh, 100));
	TEST_ASSERT(0 == fclose(fh));

	/* Do N-gram expansion on it. */
	ms_lattice_expand(l, lm);

	/* Run forward-backward on it. */
	ms_lattice_forward(l, 10);
	ms_lattice_backward(l, 10);

	TEST_ASSERT(fh = fopen("expanded.dot", "w"));
	TEST_ASSERT(0 == ms_lattice_write_dot(l, fh));
	TEST_ASSERT(0 == fclose(fh));

	/* Verify that posteriors sum to one. */

	ms_lattice_free(l);
	ngram_model_free(lm);
	logmath_free(lmath);
	return 0;
}
