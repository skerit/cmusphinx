#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <multisphinx/ms_lattice.h>

#include "test_macros.h"

int
main(int argc, char *argv[])
{
	ms_lattice_t *l;
	logmath_t *lmath;
	FILE *fh;

	lmath = logmath_init(1.0001, 0, FALSE);
	TEST_ASSERT(l = ms_lattice_init(lmath, NULL));

	/* Read in the test lattice. */

	/* Run forward on it. */

	/* Run partial backward on it for various frames. */

	/* Verify that posteriors sum to one. */

	ms_lattice_free(l);
	return 0;
}
