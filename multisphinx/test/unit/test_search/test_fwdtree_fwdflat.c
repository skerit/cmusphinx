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
	mfcc_t ***feat;
	int input_nfr, nfr, i;
	char const *hyp;
	int32 score;
	int fffr;
	char const *maxwpf = "50";

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
	acmod = ps->acmod;
	acmod2 = acmod_copy(acmod);

	cmd_ln_set_str_r(config, "-lm", TESTDATADIR "/hub4.5000.DMP");
	fwdtree = fwdtree_search_init(config, acmod, ps->dict, ps->d2p);
	fwdflat = fwdflat_search_init(config, acmod2, ps->dict, ps->d2p,
				      ((fwdtree_search_t *)fwdtree)->bptbl);

	input_nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
			      ".mfc", 0, -1, NULL, -1);
	feat = feat_array_alloc(acmod->fcb, input_nfr);
	if ((input_nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
				   ".mfc", 0, -1, feat, -1)) < 0)
		E_FATAL("Failed to read mfc file\n");
	ps_search_start(fwdtree);
	ps_search_start(fwdflat);
	i = 0;
	fffr = 0;
	for (i = 0; i < input_nfr; ++i) {
		if (fffr == input_nfr)
			fffr = 0;
		acmod_process_feat(acmod, feat[i]);
		ps_search_step(fwdtree);
		fffr += acmod_process_feat(acmod2, feat[fffr]);
		nfr = ps_search_step(fwdflat);
		E_INFO("i %d fffr %d fwdflat nfr %d\n",
		       i, fffr, nfr);
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
	feat_array_free(feat);
	ps_free(ps);

	return 0;
}
