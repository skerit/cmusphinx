#include <stdio.h>
#include <string.h>
#include <time.h>

#include <pocketsphinx.h>

#include "pocketsphinx_internal.h"
#include "acmod2.h"
#include "test_macros.h"

static int
consumer(sbthread_t *th)
{
	acmod2_t *acmod2 = sbthread_arg(th);
	int frame_idx;

	printf("Consumer %p started\n", acmod2);
	while ((frame_idx = acmod2_wait(acmod2, -1)) >= 0) {
		/* Score a frame. */
		acmod2_score(acmod2, frame_idx);
		/* Release it. */
		acmod2_release(acmod2, frame_idx);
	}
	TEST_ASSERT(acmod2_eou(acmod2));
	printf("Consumer %p exiting\n", acmod2);
	return 0;
}


int
main(int argc, char *argv[])
{
	cmd_ln_t *config;
	logmath_t *lmath;
	acmod2_t *acmod2[5];
	sbthread_t *thr[5];
	featbuf_t *fb;
	FILE *raw;
	int16 buf[2048];
	int nsamp;
	int i;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/hub4.5000.DMP",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     "-compallsen", "yes",
			     NULL);
	ps_init_defaults(config);
	fb = featbuf_init(config);
	TEST_ASSERT(fb);

	acmod2[0] = acmod2_init(config, lmath, fb);
	TEST_ASSERT(acmod2[0]);
	/* Create a couple threads to pull features out of it. */
	for (i = 0; i < 5; ++i) {
		if (i != 0)
			acmod2[i] = acmod2_copy(acmod2[0]);
		thr[i] = sbthread_start(NULL, consumer, acmod2[i]);
	}

	/* Feed them some data. */
	raw = fopen(TESTDATADIR "/chan3.raw", "rb");
	featbuf_start_utt(fb);
	while ((nsamp = fread(buf, 2, 2048, raw)) > 0) {
		int rv;
		rv = featbuf_process_raw(fb, buf, nsamp, FALSE);
		printf("Producer processed %d samples\n", nsamp);
		TEST_ASSERT(rv == 0);
	}
	fclose(raw);
	featbuf_end_utt(fb, -1);

	/* Reap those threads. */
	for (i = 0; i < 5; ++i) {
		sbthread_wait(thr[i]);
		sbthread_free(thr[i]);
	}
	cmd_ln_free_r(config);
	return 0;
}
