/**
 * @file test_garray.c Test arrays
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include "garray.h"
#include "test_macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int
main(int argc, char *argv[])
{
	garray_t *gar, *gar2;
	size_t n;
	int i;

	gar = garray_init(10, sizeof(int));
	for (i = 0; i < 10; ++i)
		garray_ent(gar, int, i) = i + 42;
	for (i = 0; i < 10; ++i)
		TEST_ASSERT(garray_ent(gar, int, i) == i + 42);
	gar2 = garray_retain(gar);
	garray_free(gar);

	n = garray_shift(gar2, 3);
	TEST_ASSERT(n == 7);
	TEST_ASSERT(garray_size(gar2) == 7);
	for (i = 0; i < 7; ++i)
		TEST_ASSERT(garray_ent(gar2, int, i) == i + 45);
	n = garray_pop(gar2, 3);
	TEST_ASSERT(n == 4);
	TEST_ASSERT(garray_size(gar2) == 4);
	for (i = 0; i < 4; ++i)
		TEST_ASSERT(garray_ent(gar2, int, i) == i + 45);
	i = 99;
	garray_append(gar, &i);
	TEST_ASSERT(garray_size(gar2) == 5);
	TEST_ASSERT(garray_ent(gar2, int, 4) == 99);
	for (i = 0; i < 10; ++i)
		garray_append(gar2, &i);
	TEST_ASSERT(garray_size(gar2) == 15);
	TEST_ASSERT(garray_ent(gar2, int, 14) == 9);
	garray_free(gar2);
	return 0;
}
