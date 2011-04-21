#include <stdio.h>
#include <string.h>
#include <time.h>

#include <multisphinx/acmod.h>
#include <multisphinx/pocketsphinx_internal.h>

#include "test_macros.h"

static int
consumer(sbthread_t *th)
{
	acmod_t *acmod = sbthread_arg(th);
	int frame_idx;

	printf("Consumer %p started\n", acmod);
	acmod_consumer_start_utt(acmod, -1);
	while ((frame_idx = acmod_consumer_wait(acmod, -1)) >= 0) {
		int senid, score;
		/* Score a frame. */
		acmod_score(acmod, frame_idx);
		score = acmod_best_score(acmod, &senid);
		printf("Consumer %p scored frame %d best %d score %d\n",
		       acmod, frame_idx, senid, score);
		/* Release it. */
		acmod_consumer_release(acmod, frame_idx);
	}
	TEST_ASSERT(acmod_eou(acmod));
	acmod_consumer_end_utt(acmod);
	printf("Consumer %p exiting\n", acmod);
	return 0;
}

int
main(int argc, char *argv[])
{
	cmd_ln_t *config;
	logmath_t *lmath;
	acmod_t *acmod[5];
	sbthread_t *thr[5];
	featbuf_t *fb;
	FILE *raw;
	int16 buf[2048];
	int nsamp;
	int i;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/bn10000.3g.arpa",
			     "-dict", TESTDATADIR "/bn10000.dic",
			     "-compallsen", "yes",
			     NULL);
	ps_init_defaults(config);
	fb = featbuf_init(config);
	TEST_ASSERT(fb);

	lmath = logmath_init(cmd_ln_float32_r(config, "-logbase"),
			     0, FALSE);
	acmod[0] = acmod_init(config, lmath, fb);
	TEST_ASSERT(acmod[0]);
	/* Create a couple threads to pull features out of it. */
	for (i = 0; i < 5; ++i) {
		if (i != 0)
			acmod[i] = acmod_copy(acmod[0]);
		thr[i] = sbthread_start(NULL, consumer, acmod[i]);
	}

	/* Feed them some data. */
	raw = fopen(TESTDATADIR "/chan3.raw", "rb");
	featbuf_producer_start_utt(fb, "chan3");
	while ((nsamp = fread(buf, 2, 2048, raw)) > 0) {
		int rv;
		rv = featbuf_producer_process_raw(fb, buf, nsamp, FALSE);
		printf("Producer processed %d samples\n", nsamp);
		TEST_ASSERT(rv > 0);
	}
	fclose(raw);
	printf("Waiting for consumers\n");
	featbuf_producer_end_utt(fb);
	printf("Finished waiting\n");

	/* Reap those threads. */
	for (i = 0; i < 5; ++i) {
		sbthread_wait(thr[i]);
		sbthread_free(thr[i]);
		acmod_free(acmod[i]);
		printf("Reaped consumer %p\n", acmod[i]);
	}
	featbuf_free(fb);
	logmath_free(lmath);
	cmd_ln_free_r(config);
	return 0;
}
