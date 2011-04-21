#include <stdio.h>

#include <sphinxbase/gq.h>
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	gq_t *q;
	int i;

	q = gq_init(sizeof(int));
	for (i = 0; i < 36; ++i)
		gq_prepend(q, &i);
	printf("%d %d\n", gq_head(q, int), gq_tail(q, int));
	TEST_ASSERT(gq_head(q, int) == 35);
	for (i = 0; i < 36; ++i) {
		gq_append(q, &i);
		printf("%d %d\n", gq_head(q, int), gq_tail(q, int));
		TEST_ASSERT(gq_head(q, int) == 35);
	}
	TEST_ASSERT(gq_tail(q, int) == 35);
	TEST_ASSERT(gq_size(q) == 72);
	for (i = 35; i >= 0; --i) {
		printf("%d %d\n", gq_head(q, int), gq_tail(q, int));
		TEST_ASSERT(gq_tail(q, int) == i);
		gq_pop(q, 1);
		TEST_ASSERT(gq_head(q, int) == i);
		gq_shift(q, 1);
	}
	gq_free(q);
	return 0;
}
