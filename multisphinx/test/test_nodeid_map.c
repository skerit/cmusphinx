#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <multisphinx/nodeid_map.h>

#include "test_macros.h"

int
main(int argc, char *argv[])
{
	nodeid_map_t *nmap;
	int fr, idx;

	srand(time(NULL)); /* Not random! */
	nmap = nodeid_map_init();
	idx = 0;
	for (fr = 0; fr < 100; ++fr) {
		int narc = rand() % 25;
		int lmstate = rand() % 100;
		int i;
		printf("%d:", fr);
		for (i = 0; i < narc; ++i) {
			printf(" %d=%d", lmstate, idx);
			nodeid_map_add(nmap, fr, lmstate++, idx++);
		}
		printf("\n");
	}
	for (fr = 0; fr < 100; ++fr) {
		nodeid_iter_t *itor;
		printf("%d:", fr);
		for (itor = nodeid_map_iter(nmap, fr);
		     itor; itor = nodeid_iter_next(itor)) {
			int32 idx, lmstate, idx2;
			idx = nodeid_iter_get(itor, &lmstate);
			idx2 = nodeid_map_map(nmap, fr, lmstate);
			printf(" %d=%d/%d", lmstate, idx, idx2);
			TEST_ASSERT(idx == idx2);
		}
		printf("\n");
	}
	nodeid_map_free(nmap);

	return 0;
}
