#include <stdio.h>
#include <string.h>
#include <time.h>

#include "pocketsphinx_internal.h"
#include "bptbl.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	ps_decoder_t *ps;
	cmd_ln_t *config;
	bptbl_t *bptbl;
	bp_t *bp;
	int fi;

	/* Get the API to initialize a bunch of stuff for us (but not the search). */
	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     "-fwdtree", "no",
			     "-fwdflat", "no",
			     "-bestpath", "no", NULL);
	ps = ps_init(config);
	bptbl = bptbl_init(ps->d2p, 10, 10);

	/* Enter a few bps starting at frame zero. */
	fi = bptbl_push_frame(bptbl, NO_BP);
	TEST_ASSERT(fi == 0);
	bp = bptbl_enter(bptbl, 42, NO_BP, 1, 0);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 0);
	TEST_ASSERT(bptbl_sf(bptbl, 0) == 0);

	fi = bptbl_push_frame(bptbl, NO_BP);
	TEST_ASSERT(fi == 1);
	bp = bptbl_enter(bptbl, 42, NO_BP, 2, 0);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 1);
	TEST_ASSERT(bptbl_sf(bptbl, 1) == 0);

	fi = bptbl_push_frame(bptbl, NO_BP);
	TEST_ASSERT(fi == 2);
	bp = bptbl_enter(bptbl, 42, NO_BP, 3, 0);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 2);
	TEST_ASSERT(bptbl_sf(bptbl, 2) == 0);

	fi = bptbl_push_frame(bptbl, NO_BP);
	TEST_ASSERT(fi == 3);
	bp = bptbl_enter(bptbl, 69, 1, 4, 0);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 3);
	TEST_ASSERT(bptbl_sf(bptbl, 3) == 2);
	bp = bptbl_enter(bptbl, 69, 1, 5, 0);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 4);
	TEST_ASSERT(bptbl_sf(bptbl, 4) == 2);

	dump_bptable(bptbl);
	/* This should cause frames 0 and 1 to get garbage collected,
	 * invalidating bp #0 and renumbering bp #1 to 0.  Ensure that
	 * everything else is still the same.
	 */
	fi = bptbl_push_frame(bptbl, 2);
	TEST_ASSERT(fi == 4);
	dump_bptable(bptbl);

	/* This one is retired. */
	bp = bptbl_ent(bptbl, 0);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 0);
	TEST_ASSERT(bptbl_sf(bptbl, 0) == 0);
	TEST_ASSERT(bp->wid == 42);
	TEST_ASSERT(bp->score == 2);

	/* FIXME: bptbl_ent(bptbl, 1) should return NULL since it is
	 * now an invalid index. */

	/* This one is the first active one.  It has not been renumbered. */
	bp = bptbl_ent(bptbl, 2);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 2);
	TEST_ASSERT(bptbl_sf(bptbl, 2) == 0);
	TEST_ASSERT(bp->wid == 42);
	TEST_ASSERT(bp->score == 3);

	bp = bptbl_ent(bptbl, 3);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 3);
	TEST_ASSERT(bptbl_sf(bptbl, 3) == 2);
	TEST_ASSERT(bp->wid == 69);
	TEST_ASSERT(bp->score == 4);

	bp = bptbl_ent(bptbl, 4);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 4);
	TEST_ASSERT(bptbl_sf(bptbl, 4) == 2);
	TEST_ASSERT(bp->wid == 69);
	TEST_ASSERT(bp->score == 5);

	/* Add some more bps and gc again. */
	fi = bptbl_push_frame(bptbl, 2);
	TEST_ASSERT(fi == 5);
	bp = bptbl_enter(bptbl, 999, 3, 5, 0);
	TEST_ASSERT(bptbl_idx(bptbl, bp) == 5);
	TEST_ASSERT(bptbl_sf(bptbl, 5) == 4);

	dump_bptable(bptbl);
	/* This should cause frames 2 through 4 to get garbage
	 * collected.
	 */
	fi = bptbl_push_frame(bptbl, 5);
	TEST_ASSERT(fi == 6);
	dump_bptable(bptbl);

	bptbl_free(bptbl);
	ps_free(ps);

	return 0;
}
