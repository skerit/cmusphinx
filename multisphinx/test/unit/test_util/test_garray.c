/**
 * @file test_garray.c Test arrays
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include "garray.h"
#include "test_macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int
test_indexing(void)
{
	garray_t *gar, *gar2;
	size_t n;
	int i;

	gar = garray_init(10, sizeof(int));
	for (i = 0; i < 10; ++i)
		garray_ent(gar, int, i) = i + 42;
	for (i = 0; i < 10; ++i)
		TEST_ASSERT(garray_ent(gar, int, i) == i + 42);

	garray_set_cmp(gar, garray_cmp_int32, NULL);
	i = 45;
	n = garray_find_first(gar, &i);
	TEST_ASSERT(garray_ent(gar, int, n) == 45);
	TEST_ASSERT(n == 3);

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
	gar = garray_slice(gar2, 10, 5);
	garray_free(gar2);
	TEST_ASSERT(garray_ent(gar, int, 4) == 9);

	TEST_ASSERT(garray_set_base(gar, 68) == 0);
	TEST_ASSERT(garray_ent(gar, int, 72) == 9);

	garray_free(gar);
	return 0;
}

static int
test_sorting(void)
{
	static char const *data[] = {
		"eggs",
		"spam",
		"bacon",
		"eggs",
		"spam",
		"spam",
		"SPAM",
		"potatoes",
		"pie"
	};
	garray_t *gar;
	int i;

	gar = garray_init(0, sizeof(char const *));
	for (i = 0; i < 9; ++i)
		garray_append(gar, &data[i]);

	garray_set_cmp(gar, garray_cmp_str, NULL);
	garray_sort(gar);

	for (i = 0; i < 8; ++i) {
		printf("%s\n", garray_ent(gar, char const *, i));
		TEST_ASSERT(strcmp(garray_ent(gar, char const *, i),
				   garray_ent(gar, char const *, i)) <= 0);
	}
	printf("%s\n", garray_ent(gar, char const *, 8));
	garray_free(gar);
	return 0;
}

static int
test_insertion(void)
{
	garray_t *gar;
	int i;

	gar = garray_init(10, sizeof(int));
	for (i = 0; i < 10; ++i)
		garray_ent(gar, int, i) = i + 1;
	i = 0;
	garray_insert(gar, 0, &i);
	for (i = 0; i <= 10; ++i)
		TEST_ASSERT(i == garray_ent(gar, int, i));
	garray_free(gar);
	return 0;
}

static int
test_deletion(void)
{
	garray_t *gar;
	int i;

	gar = garray_init(15, sizeof(int));
	for (i = 0; i < 15; ++i)
		garray_ent(gar, int, i) = i;
	garray_delete(gar, 0, 5);
	for (i = 0; i < 10; ++i)
		TEST_ASSERT(i + 5 == garray_ent(gar, int, i));
	garray_free(gar);
	return 0;
}

int
main(int argc, char *argv[])
{
	test_indexing();
	test_sorting();
	test_insertion();
	test_deletion();
	return 0;
}
