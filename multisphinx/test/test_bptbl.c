#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/strfuncs.h>
#include <sphinxbase/sbthread.h>

#include <multisphinx/bptbl.h>
#include <multisphinx/search.h>

#include "test_macros.h"

int
main(int argc, char *argv[])
{
	bin_mdef_t *mdef;
	dict2pid_t *d2p;
	dict_t *dict;
	cmd_ln_t *config;
	seg_iter_t *seg;
	bptbl_t *bptbl;
	bp_t bpe;
	bpidx_t bp;
	int fi;
	int i;

	/* Get the API to initialize a bunch of stuff for us (but not the search). */
	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/bn10000.3g.arpa",
			     "-dict", TESTDATADIR "/bn10000.dic",
			     NULL);
	ps_init_defaults(config);
	mdef = bin_mdef_read(config, cmd_ln_str_r(config, "-mdef"));
	dict = dict_init(config, mdef);
	d2p = dict2pid_build(mdef, dict);

	bptbl = bptbl_init("test", d2p, 10, 10);

	/* Enter a few bps starting at frame zero. */
	fi = bptbl_push_frame(bptbl, NO_BP);
	TEST_ASSERT(fi == 0);
	bp = bptbl_enter(bptbl, 42, NO_BP, 1, 0);
	TEST_ASSERT(bp == 0);
	TEST_ASSERT(bptbl_sf(bptbl, bp) == 0);

	fi = bptbl_push_frame(bptbl, NO_BP);
	TEST_ASSERT(fi == 1);
	bp = bptbl_enter(bptbl, 42, NO_BP, 2, 0);
	TEST_ASSERT(bp == 1);
	TEST_ASSERT(bptbl_sf(bptbl, bp) == 0);

	fi = bptbl_push_frame(bptbl, NO_BP);
	TEST_ASSERT(fi == 2);
	bp = bptbl_enter(bptbl, 42, NO_BP, 3, 0);
	TEST_ASSERT(bp == 2);
	TEST_ASSERT(bptbl_sf(bptbl, bp) == 0);

	fi = bptbl_push_frame(bptbl, NO_BP);
	TEST_ASSERT(fi == 3);
	bp = bptbl_enter(bptbl, 69, 1, 4, 0);
	TEST_ASSERT(bp == 3);
	TEST_ASSERT(bptbl_sf(bptbl, bp) == 2);
	bp = bptbl_enter(bptbl, 69, 1, 5, 0);
	TEST_ASSERT(bp == 4);
	TEST_ASSERT(bptbl_sf(bptbl, bp) == 2);

	bptbl_dump(bptbl);
	/* This should cause frames 0 and 1 to get garbage collected,
	 * invalidating bp #0 and renumbering bp #1 to 0.  Ensure that
	 * everything else is still the same.
	 */
	fi = bptbl_push_frame(bptbl, 2);
	TEST_ASSERT(fi == 4);
	bptbl_dump(bptbl);

	/* This one is retired. */
	bptbl_get_bp(bptbl, 0, &bpe);
	TEST_ASSERT(bptbl_sf(bptbl, 0) == 0);
	TEST_ASSERT(bpe.wid == 42);
	TEST_ASSERT(bpe.score == 2);

	/* FIXME: bptbl_ent(bptbl, 1) should return NULL since it is
	 * now an invalid index. */

	/* This one is the first active one.  It has not been renumbered. */
	bptbl_get_bp(bptbl, 2, &bpe);
	TEST_ASSERT(bptbl_sf(bptbl, 2) == 0);
	TEST_ASSERT(bpe.wid == 42);
	TEST_ASSERT(bpe.score == 3);

	bptbl_get_bp(bptbl, 3, &bpe);
	TEST_ASSERT(bptbl_sf(bptbl, 3) == 2);
	TEST_ASSERT(bpe.wid == 69);
	TEST_ASSERT(bpe.score == 4);

	bptbl_get_bp(bptbl, 4, &bpe);
	TEST_ASSERT(bptbl_sf(bptbl, 4) == 2);
	TEST_ASSERT(bpe.wid == 69);
	TEST_ASSERT(bpe.score == 5);

	/* Add some more bps and gc again. */
	fi = bptbl_push_frame(bptbl, 2);
	TEST_ASSERT(fi == 5);
	bp = bptbl_enter(bptbl, 999, 3, 5, 0);
	TEST_ASSERT(bp == 5);
	TEST_ASSERT(bptbl_sf(bptbl, bp) == 4);

	bptbl_dump(bptbl);
	/* This should cause frames 2 through 4 to get garbage
	 * collected.
	 */
	fi = bptbl_push_frame(bptbl, 5);
	TEST_ASSERT(fi == 6);
	bptbl_dump(bptbl);

	/* Now add a bunch of stuff to see what happens. */
	for (i = 0; i < 6; ++i) {
		bp = bptbl_enter(bptbl, 42, 5, 6 + i, 0);
	}
	fi = bptbl_push_frame(bptbl, 6);
	TEST_ASSERT(fi == 7);
	for (i = 0; i < 3; ++i) {
		bp = bptbl_enter(bptbl, 69, 6, 12 + i, 0);
	}

	/* Finalize the backpointer table (i.e. retire all active) */
	bptbl_dump(bptbl);
	bptbl_finalize(bptbl);
	bptbl_dump(bptbl);

	/* FIXME: assert the correct values here. */
	for (seg = bptbl_seg_iter(bptbl, NULL, BAD_S3WID); seg;
	     seg = ps_seg_next(seg)) {
		char const *word;
		int sf, ef;
		int32 post, lscr, ascr, lback;

		word = ps_seg_word(seg);
		ps_seg_frames(seg, &sf, &ef);
		post = ps_seg_prob(seg, &ascr, &lscr, &lback);
		printf("%s (%d:%d) ascr = %d lscr = %d lback = %d\n",
		       word, sf, ef,
		       ascr, lscr, lback);
	}

	bptbl_free(bptbl);
	dict2pid_free(d2p);
	dict_free(dict);
	bin_mdef_free(mdef);

	return 0;
}
