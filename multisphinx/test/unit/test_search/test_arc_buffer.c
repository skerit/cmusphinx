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
	int fi, i, next_sf;
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

	arcs = fwdflat_arc_buffer_init();

	/* Enter a bunch of initial bps (like silence) */
	fi = bptbl_push_frame(bptbl, NO_BP);
	bp = bptbl_enter(bptbl, 42, NO_BP, 1, 0);
	fi = bptbl_push_frame(bptbl, NO_BP);
	bp = bptbl_enter(bptbl, 42, NO_BP, 2, 0);
	fi = bptbl_push_frame(bptbl, NO_BP);
	bp = bptbl_enter(bptbl, 42, NO_BP, 3, 0);

	/* Enter a couple of words pointing back to the silences. */
	fi = bptbl_push_frame(bptbl, NO_BP);
	bp = bptbl_enter(bptbl, 69, 1, 4, 0);
	bp = bptbl_enter(bptbl, 69, 1, 5, 0);

	/* Garbage collect some things. */
	fi = bptbl_push_frame(bptbl, 2);

	/* Add some more words. */
	fi = bptbl_push_frame(bptbl, 2);
	bp = bptbl_enter(bptbl, 999, 4, 5, 0);

	/* Garbage collect some things. */
	fi = bptbl_push_frame(bptbl, 5);
	bptbl_dump(bptbl);
	next_sf = bptbl_ent(bptbl,
			    bptbl->oldest_bp)->frame + 1;
	E_INFO("next_sf %d\n", next_sf);
	fwdflat_arc_buffer_extend(arcs, next_sf);
	i = fwdflat_arc_buffer_add_bps(arcs, bptbl,
				       0, bptbl->first_invert_bp);
	E_INFO("Added %d arcs\n", i);
	fwdflat_arc_buffer_commit(arcs);

	/* Now add a bunch of stuff to see what happens. */
	for (i = 0; i < 6; ++i) {
		bp = bptbl_enter(bptbl, 42, 5, 6 + i, 0);
	}
	fi = bptbl_push_frame(bptbl, 9);
	for (i = 0; i < 3; ++i) {
		bp = bptbl_enter(bptbl, 69, 6, 12 + i, 0);
	}
	fi = bptbl_push_frame(bptbl, 12);
	bptbl_dump(bptbl);
	next_sf = bptbl_ent(bptbl,
			    bptbl->oldest_bp)->frame + 1;
	E_INFO("next_sf %d\n", next_sf);
	fwdflat_arc_buffer_extend(arcs, next_sf);
	i = fwdflat_arc_buffer_add_bps(arcs, bptbl,
				       0, bptbl->first_invert_bp);
	E_INFO("Added %d arcs\n", i);
	fwdflat_arc_buffer_commit(arcs);

	for (i = 0; i < 3; ++i) {
		bp = bptbl_enter(bptbl, 420, 6, 39 + i, 0);
	}
	bptbl_finalize(bptbl);
	bptbl_dump(bptbl);
	next_sf = bptbl_ent(bptbl,
			    bptbl->oldest_bp)->frame + 1;
	E_INFO("next_sf %d\n", next_sf);
	fwdflat_arc_buffer_extend(arcs, next_sf);
	i = fwdflat_arc_buffer_add_bps(arcs, bptbl,
				       0, bptbl->first_invert_bp);
	E_INFO("Added %d arcs\n", i);
	fwdflat_arc_buffer_commit(arcs);

	fwdflat_arc_buffer_free(arcs);
	bptbl_free(bptbl);
	ps_free(ps);
	return 0;
}
