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
	bin_mdef_t *mdef;
	dict2pid_t *d2p;
	dict_t *dict;
	logmath_t *lmath;
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
			     "-lm", TESTDATADIR "/hub4.5000.DMP",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     "-fwdflatefwid", "3",
			     "-maxwpf", maxwpf,
			     "-latsize", "512", NULL);

	ps_init_defaults(config);
	lmath = logmath_init(cmd_ln_float32_r(config, "-logbase"),
			     0, FALSE);

	acmod = acmod_init(config, lmath, NULL, NULL);
	acmod2 = acmod_copy(acmod);
	mdef = bin_mdef_read(config, cmd_ln_str_r(config, "-mdef"));
	dict = dict_init(config, mdef);
	d2p = dict2pid_build(mdef, dict);

	fwdtree = fwdtree_search_init(config, acmod, dict, d2p);
	fwdflat = fwdflat_search_init(config, acmod2, dict, d2p,
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
#if 0 /* FIXME: Doesn't work yet.  Splice results together. */
	hyp = ps_search_hyp(fwdflat, &score);
	printf("fwdflat: %s (%d)\n", hyp, score);
#endif

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

	acmod_free(acmod);
	acmod_free(acmod2);
	ps_search_free(fwdtree);
	ps_search_free(fwdflat);
	feat_array_free(feat);
	logmath_free(lmath);
	cmd_ln_free_r(config);

	return 0;
}
