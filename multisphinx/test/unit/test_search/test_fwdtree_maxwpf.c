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
	acmod_t *acmod;
	ps_search_t *fwdtree;
	mfcc_t ***feat;
	int nfr, i, k;
	char const *hyp;
	int32 score;
	char const *maxwpf = "50";
	int ncat = 5;

	if (argc > 1)
		maxwpf = argv[1];
	if (argc > 2)
		ncat = atoi(argv[2]);
	ckd_set_jump(NULL, TRUE);
	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     "-fwdtree", "no",
			     "-fwdflat", "no",
			     "-maxwpf", maxwpf,
			     "-latsize", "512",
			     "-bestpath", "no", NULL);

	/* Get the API to initialize a bunch of stuff for us (but not the search). */
	ps = ps_init(config);
	acmod = ps->acmod;
	cmd_ln_set_str_r(config, "-lm", TESTDATADIR "/hub4.5000.DMP");
	fwdtree = fwdtree_search_init(config, acmod, ps->dict, ps->d2p);

	nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
			      ".mfc", 0, -1, NULL, -1);
	feat = feat_array_alloc(acmod->fcb, nfr);
	if ((nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
				   ".mfc", 0, -1, feat, -1)) < 0)
		E_FATAL("Failed to read mfc file\n");
	ps_search_start(fwdtree);
	for (k = 0; k < ncat; ++k) {
		for (i = 0; i < nfr; ++i) {
			acmod_process_feat(acmod, feat[i]);
			ps_search_step(fwdtree);
		}
	}
	ps_search_finish(fwdtree);
	hyp = ps_search_hyp(fwdtree, &score);
	printf("hyp: %s (%d)\n", hyp, score);

	ps_search_free(fwdtree);
	feat_array_free(feat);
	ps_free(ps);

	return 0;
}
