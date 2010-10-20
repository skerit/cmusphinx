#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/sbthread.h>
#include <pocketsphinx.h>

#include "pocketsphinx_internal.h"
#include "test_macros.h"

/* This will be a prototype for how acmod works. */
typedef struct feat_reader_s {
	featbuf_t *src;
	int done;
	int fr;
} feat_reader_t;

static int
get_some_features(sbthread_t *th)
{
	feat_reader_t *fr = sbthread_arg(th);

	while (!fr->done) {
		/* Wait for frame to be available. */
		featbuf_wait(fr->src, fr->fr, 1000);

		/* Do something with that frame. */

		/* Release that frame. */
		featbuf_release(fr->src, fr->fr);
		++fr->fr;
	}
	return 0;
}

int
main(int argc, char *argv[])
{
	cmd_ln_t *config;
	featbuf_t *fb;
	feat_t *feat;
	fe_t *fe;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/hub4.5000.DMP",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     NULL);
	ps_init_defaults(config);
	fb = featbuf_init(config);
	TEST_ASSERT(fb);

	/* Assert that the feature parameters are correct. */
	fe = featbuf_get_fe(fb);
	feat = featbuf_get_fcb(fb);

	/* Create a couple threads to pull features out of it. */

	/* Feed it some data. */

	/* Reap those threads. */
	featbuf_free(fb);
	cmd_ln_free(config);
	return 0;
}
