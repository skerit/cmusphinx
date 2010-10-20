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
	bin_mdef_t *mdef;
	dict2pid_t *d2p;
	dict_t *dict;
	cmd_ln_t *config;
	acmod_t *acmod, *acmod2;
	ps_search_t *fwdtree, *fwdtree2;
	mfcc_t ***feat;
	int nfr, i;
	char const *hyp, *hyp2;
	int32 score, score2;
	logmath_t *lmath;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/hub4.5000.DMP",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     NULL);
	ps_init_defaults(config);
	lmath = logmath_init(cmd_ln_float32_r(config, "-logbase"),
			     0, FALSE);

	/* FIXME: This actually only tests s2_semi_mgau for now. */
	acmod = acmod_init(config, lmath, NULL, NULL);
	acmod2 = acmod_copy(acmod);
	mdef = bin_mdef_read(config, cmd_ln_str_r(config, "-mdef"));
	dict = dict_init(config, mdef);
	d2p = dict2pid_build(mdef, dict);

	fwdtree = fwdtree_search_init(config, acmod, dict, d2p);
	fwdtree2 = fwdtree_search_init(config, acmod2, dict, d2p);

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
		ps_search_step(fwdtree);
		acmod_advance(acmod);
		ps_search_step(fwdtree2);
		acmod_advance(acmod2);
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
	acmod_free(acmod);
	acmod_free(acmod2);
	logmath_free(lmath);
	cmd_ln_free_r(config);

	return 0;
}
