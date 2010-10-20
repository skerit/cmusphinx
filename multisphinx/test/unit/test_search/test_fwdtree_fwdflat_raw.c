#include <pocketsphinx.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/feat.h>
#include <sphinxbase/ckd_alloc.h>

#include "pocketsphinx_internal.h"
#include "fwdtree_search.h"
#include "fwdflat_search.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	ps_decoder_t *ps;
	cmd_ln_t *config;
	acmod_t *acmod, *acmod2;
	ps_search_t *fwdtree;
	ps_search_t *fwdflat;
	int nfr, i;
	char const *hyp;
	int32 score;
	char const *maxwpf = "50";
	FILE *rawfh, *rawfh2;
	int16 buf[2048], buf2[2048];
	int16 const *bptr, *bptr2;

	if (argc > 1)
		maxwpf = argv[1];
	ckd_set_jump(NULL, TRUE);
	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     "-fwdtree", "no",
			     "-fwdflat", "no",
			     "-fwdflatefwid", "3",
			     "-maxwpf", maxwpf,
			     "-latsize", "512",
			     "-bestpath", "no", NULL);

	/* Get the API to initialize a bunch of stuff for us (but not the search). */
	ps = ps_init(config);

	/* Create two acmods - this is the first take at running stuff
	 * from a file in parallel. */
	acmod = ps->acmod;
	acmod2 = acmod_copy(acmod);

	cmd_ln_set_str_r(config, "-lm", TESTDATADIR "/hub4.5000.DMP");
	fwdtree = fwdtree_search_init(config, acmod, ps->dict, ps->d2p);
	fwdflat = fwdflat_search_init(config, acmod2, ps->dict, ps->d2p,
				      ((fwdtree_search_t *)fwdtree)->bptbl);

	TEST_ASSERT(rawfh = fopen(TESTDATADIR "/goforward.raw", "rb"));
	TEST_ASSERT(rawfh2 = fopen(TESTDATADIR "/goforward.raw", "rb"));
	/* FIXME: I think the search module should probably do this itself. */
	TEST_EQUAL(0, acmod_start_utt(acmod));
	TEST_EQUAL(0, acmod_start_utt(acmod2));
	ps_search_start(fwdtree);
	ps_search_start(fwdflat);
	while (!feof(rawfh)) {
		nread = fread(buf, sizeof(*buf), 2048, rawfh);
		bptr = buf;
		while ((nfr = acmod_process_raw(acmod, &bptr, &nread, FALSE)) > 0) {
			while (acmod->n_feat_frame > 0) {
				ps_search_step(fwdtree);
			}
		}
		nread = fread(buf2, sizeof(*buf2), 2048, rawfh2);
		bptr2 = buf2;
		while ((nfr = acmod_process_raw(acmod, &bptr2, &nread, FALSE)) > 0) {
			while (acmod->n_feat_frame > 0) {
				ps_search_step(fwdflat);
			}
		}
	}
	ps_search_finish(fwdtree);
	hyp = ps_search_hyp(fwdtree, &score);
	printf("fwdtree: %s (%d)\n", hyp, score);
	while (fffr < input_nfr) {
		fffr += acmod_process_feat(acmod2, feat[fffr]);
		if ((nfr = ps_search_step(fwdflat)) <= 0)
			break;
		E_INFO("i %d fffr %d fwdflat nfr %d\n",
		       i, fffr, nfr);
	}
	ps_search_finish(fwdflat);
	hyp = ps_search_hyp(fwdflat, &score);
	printf("fwdflat: %s (%d)\n", hyp, score);

	acmod_free(acmod2);
	ps_search_free(fwdtree);
	ps_search_free(fwdflat);
	ps_free(ps);

	return 0;
}
