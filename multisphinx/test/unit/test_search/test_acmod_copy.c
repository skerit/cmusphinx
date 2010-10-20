#include <pocketsphinx.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/feat.h>

#include "pocketsphinx_internal.h"
#include "fwdtree_search.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	ps_decoder_t *ps;
	cmd_ln_t *config;
	acmod_t *acmod, *acmod2;
	ps_search_t *fwdtree, *fwdtree2;
	mfcc_t ***feat;
	int nfr, i;
	char const *hyp, *hyp2;
	int32 score, score2;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     "-fwdtree", "no",
			     "-fwdflat", "no",
			     "-bestpath", "no", NULL);

	/* Get the API to initialize a bunch of stuff for us (but not the search). */
	ps = ps_init(config);
	/* FIXME: This actually only tests s2_semi_mgau for now. */
	acmod = ps->acmod;
	acmod2 = acmod_copy(acmod);

	cmd_ln_set_str_r(config, "-lm", TESTDATADIR "/hub4.5000.DMP");
	fwdtree = fwdtree_search_init(config, acmod, ps->dict, ps->d2p);
	fwdtree2 = fwdtree_search_init(config, acmod2, ps->dict, ps->d2p);

	nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
			      ".mfc", 0, -1, NULL, -1);
	feat = feat_array_alloc(acmod->fcb, nfr);
	if ((nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
				   ".mfc", 0, -1, feat, -1)) < 0)
		E_FATAL("Failed to read mfc file\n");

	ps_search_start(fwdtree);
	ps_search_start(fwdtree2);
	for (i = 0; i < nfr; ++i) {
		acmod_process_feat(acmod, feat[i]);
		acmod_process_feat(acmod2, feat[i]);
		while (acmod->n_feat_frame > 0) {
			ps_search_step(fwdtree, acmod->output_frame);
			acmod_advance(acmod);
		}
		while (acmod2->n_feat_frame > 0) {
			ps_search_step(fwdtree2, acmod2->output_frame);
			acmod_advance(acmod2);
		}
	}
	ps_search_finish(fwdtree);
	hyp = ps_search_hyp(fwdtree, &score);
	printf("hyp: %s (%d)\n", hyp, score);

	ps_search_finish(fwdtree2);
	hyp2 = ps_search_hyp(fwdtree2, &score2);
	printf("hyp2: %s (%d)\n", hyp2, score2);
	TEST_ASSERT(score == score2);
	TEST_ASSERT(0 == strcmp(hyp, hyp2));

	ps_search_free(fwdtree);
	ps_search_free(fwdtree2);
	feat_array_free(feat);
	acmod_free(acmod2);
	ps_free(ps);

	return 0;
}
