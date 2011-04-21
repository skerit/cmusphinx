#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/sbthread.h>

#include <multisphinx/pocketsphinx_internal.h>

#include "test_macros.h"

/* This will be a prototype for how acmod works. */
typedef struct feat_reader_s {
	featbuf_t *src;
	int fr;
} feat_reader_t;

static int
consumer(sbthread_t *th)
{
	feat_reader_t *fr = sbthread_arg(th);
	mfcc_t feat[52];

	printf("Consumer %p started\n", fr);
	while (1) {
		/* Wait for frame to be available. */
		if (featbuf_consumer_wait(fr->src, fr->fr, -1, feat) < 0)
			break;

		/* Do something with that frame. */
		printf("Consumer %p got frame %d: %f %f\n",
		       fr, fr->fr, feat[0], feat[1]);

		/* Release that frame. */
		featbuf_consumer_release(fr->src, fr->fr, fr->fr + 1);
		++fr->fr;
	}
	featbuf_consumer_end_utt(fr->src, fr->fr);
	printf("Consumer %p done\n", fr);
	sleep(2);
	printf("Consumer %p exiting\n", fr);
	return 0;
}

int
main(int argc, char *argv[])
{
	cmd_ln_t *config;
	feat_reader_t *fr[5];
	sbthread_t *thr[5];
	featbuf_t *fb;
	feat_t *feat;
	fe_t *fe;
	FILE *raw;
	int16 buf[2048];
	int nsamp;
	int i;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/bn10000.3g.arpa",
			     "-dict", TESTDATADIR "/bn10000.dic",
			     NULL);
	ps_init_defaults(config);
	fb = featbuf_init(config);
	TEST_ASSERT(fb);

	/* Assert that the feature parameters are correct. */
	fe = featbuf_get_fe(fb);
	feat = featbuf_get_fcb(fb);

	/* Create a couple threads to pull features out of it. */
	for (i = 0; i < 5; ++i) {
		fr[i] = ckd_calloc(1, sizeof(**fr));
		fr[i]->src = featbuf_retain(fb);
		thr[i] = sbthread_start(NULL, consumer, fr[i]);
	}

	/* Feed them some data. */
	raw = fopen(TESTDATADIR "/chan3.raw", "rb");
	featbuf_producer_start_utt(fb, "utt");
	while ((nsamp = fread(buf, 2, 2048, raw)) > 0) {
		int rv;
		rv = featbuf_producer_process_raw(fb, buf, nsamp, FALSE);
		printf("Producer processed %d samples: %d\n", nsamp, rv);
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
		printf("featbuf rc %d\n", featbuf_free(fr[i]->src));
		ckd_free(fr[i]);
		printf("Reaped consumer %p\n", fr[i]);
	}
	printf("featbuf rc %d\n", featbuf_free(fb));
	cmd_ln_free_r(config);
	return 0;
}
