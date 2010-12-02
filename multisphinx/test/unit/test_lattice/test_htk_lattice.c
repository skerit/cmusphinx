#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "ms_lattice.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	ms_lattice_t *l;
	logmath_t *lmath;
	FILE *fh;

	lmath = logmath_init(1.0001, 0, FALSE);
	TEST_ASSERT(l = ms_lattice_init(lmath));
	TEST_ASSERT(fh = fopen(TESTDATADIR "/050c0103.slf", "r"));
	TEST_ASSERT(0 == ms_lattice_read_htk(l, fh));

	ms_lattice_free(l);
	return 0;
}
