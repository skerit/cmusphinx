#include "pocketsphinx_internal.h"
#include "fwdflat_search.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	bin_mdef_t *mdef;
	dict2pid_t *d2p;
	dict_t *dict;
	arc_buffer_t *arcs;
	arc_t *a;
	cmd_ln_t *config;
	bptbl_t *bptbl;
	int fi, i, next_sf;
	bpidx_t bp;
	bp_t bpe;

	/* Get the API to initialize a bunch of stuff for us (but not the search). */
	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/hub4.5000.DMP",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     NULL);
	ps_init_defaults(config);
	mdef = bin_mdef_read(config, cmd_ln_str_r(config, "-mdef"));
	dict = dict_init(config, mdef);
	d2p = dict2pid_build(mdef, dict);

	bptbl = bptbl_init(d2p, 10, 10);

	arcs = arc_buffer_init(NULL);

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
	next_sf = bptbl_active_sf(bptbl);
	E_INFO("next_sf %d\n", next_sf);
	arc_buffer_extend(arcs, next_sf);
	i = arc_buffer_add_bps(arcs, bptbl,
			       0, bptbl_retired_idx(bptbl));
	E_INFO("Added %d arcs\n", i);
	arc_buffer_commit(arcs);

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
	bptbl_get_bp(bptbl, bptbl->oldest_bp, &bpe);
	next_sf = bpe.frame + 1;
	E_INFO("next_sf %d\n", next_sf);
	arc_buffer_extend(arcs, next_sf);
	i = arc_buffer_add_bps(arcs, bptbl,
				       0, bptbl_retired_idx(bptbl));
	E_INFO("Added %d arcs\n", i);
	arc_buffer_commit(arcs);

	for (i = 0; i < 3; ++i) {
		bp = bptbl_enter(bptbl, 420, 6, 39 + i, 0);
	}
	bptbl_finalize(bptbl);
	bptbl_dump(bptbl);
	bptbl_get_bp(bptbl, bptbl->oldest_bp, &bpe);
	next_sf = bpe.frame + 1;
	E_INFO("next_sf %d\n", next_sf);
	arc_buffer_extend(arcs, next_sf);
	i = arc_buffer_add_bps(arcs, bptbl,
				       0, bptbl_retired_idx(bptbl));
	E_INFO("Added %d arcs\n", i);
	arc_buffer_commit(arcs);

	a = arc_buffer_iter(arcs, 2);
	TEST_ASSERT(a->wid == 69);
	a = arc_next(arcs, a);
	TEST_ASSERT(a->wid == 69);
	a = arc_next(arcs, a);
	TEST_ASSERT(a == arc_buffer_iter(arcs, 3));

	for (a = arc_buffer_iter(arcs, 6);
	     a != arc_buffer_iter(arcs, 8);
	     a = arc_next(arcs, a)) {
		TEST_ASSERT(a->wid == 42 || a->wid == 420);
		TEST_ASSERT(a->src >= 6 && a->src < 8);
	}

	arc_buffer_dump(arcs, dict);
	arc_buffer_release(arcs, 6);
	for (a = arc_buffer_iter(arcs, 6);
	     a != arc_buffer_iter(arcs, 8);
	     a = arc_next(arcs, a)) {
		printf("%d %d %d\n", a->wid, a->src, a->dest);
		TEST_ASSERT(a->wid == 42 || a->wid == 420);
		TEST_ASSERT(a->src >= 6 && a->src < 8);
	}
	arc_buffer_dump(arcs, dict);
	arc_buffer_free(arcs);
	bptbl_free(bptbl);
	dict2pid_free(d2p);
	dict_free(dict);
	bin_mdef_free(mdef);

	return 0;
}
