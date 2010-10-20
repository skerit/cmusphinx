#include "pocketsphinx_internal.h"
#include "fwdflat_search.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	fwdflat_arc_buffer_t *arcs;
	ps_decoder_t *ps;
	cmd_ln_t *config;
	bptbl_t *bptbl;
	int fi;
	bp_t *bp;

	/* Get the API to initialize a bunch of stuff for us (but not the search). */
	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     "-fwdtree", "no",
			     "-fwdflat", "no",
			     "-bestpath", "no", NULL);
	ps = ps_init(config);
	bptbl = bptbl_init(ps->d2p, 10, 10);

	arcs = fwdflat_arc_buffer_init(3, 10);

	fi = bptbl_push_frame(bptbl, NO_BP);
	bp = bptbl_enter(bptbl, 42, NO_BP, 1, 0);
	fi = bptbl_push_frame(bptbl, NO_BP);
	bp = bptbl_enter(bptbl, 42, NO_BP, 2, 0);
	fi = bptbl_push_frame(bptbl, NO_BP);
	bp = bptbl_enter(bptbl, 42, NO_BP, 3, 0);

	fi = bptbl_push_frame(bptbl, NO_BP);
	bp = bptbl_enter(bptbl, 69, 1, 4, 0);
	bp = bptbl_enter(bptbl, 69, 1, 5, 0);
	bptbl_dump(bptbl);

	/* This should cause frames 0 and 1 to get garbage collected,
	 * invalidating bp #0 and renumbering bp #1 to 0.  Ensure that
	 * everything else is still the same.
	 */
	fi = bptbl_push_frame(bptbl, 2);
	bptbl_dump(bptbl);

	/* This one is retired. */
	bp = bptbl_ent(bptbl, 0);

	fi = bptbl_push_frame(bptbl, 2);
	bp = bptbl_enter(bptbl, 999, 3, 5, 0);

	bptbl_dump(bptbl);
	/* This should cause frames 2 through 4 to get garbage
	 * collected.
	 */
	fi = bptbl_push_frame(bptbl, 5);
	bptbl_dump(bptbl);

	/* Now add a bunch of stuff to see what happens. */
	for (i = 0; i < 6; ++i) {
		bp = bptbl_enter(bptbl, 42, 5, 6 + i, 0);
	}
	fi = bptbl_push_frame(bptbl, 6);
	for (i = 0; i < 3; ++i) {
		bp = bptbl_enter(bptbl, 69, 6, 12 + i, 0);
	}

	bptbl_free(bptbl);
	ps_free(ps);
	return 0;
}
