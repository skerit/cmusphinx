#include <pocketsphinx.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/feat.h>

#include "pocketsphinx_internal.h"
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
	acmod_t *acmod;
	ps_search_t *fwdflat;
	mfcc_t ***feat;
	int nfr, i;
	char const *hyp;
	int32 score;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/hub4.5000.DMP",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     NULL);
	ps_init_defaults(config);
	lmath = logmath_init(cmd_ln_float32_r(config, "-logbase"),
			     0, FALSE);
	acmod = acmod_init(config, lmath, NULL, NULL);
	mdef = bin_mdef_read(config, cmd_ln_str_r(config, "-mdef"));
	dict = dict_init(config, mdef);
	d2p = dict2pid_build(mdef, dict);
	fwdflat = fwdflat_search_init(config, acmod, dict, d2p, NULL);

	nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
			      ".mfc", 0, -1, NULL, -1);
	feat = feat_array_alloc(acmod->fcb, nfr);
	if ((nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
				   ".mfc", 0, -1, feat, -1)) < 0)
		E_FATAL("Failed to read mfc file\n");
	ps_search_start(fwdflat);
	for (i = 0; i < nfr; ++i) {
		acmod_process_feat(acmod, feat[i]);
		ps_search_step(fwdflat);
	}
	ps_search_finish(fwdflat);
	hyp = ps_search_hyp(fwdflat, &score);
	printf("hyp: %s (%d)\n", hyp, score);

	ps_search_free(fwdflat);
	feat_array_free(feat);
	acmod_free(acmod);
	logmath_free(lmath);
	cmd_ln_free_r(config);

	return 0;
}
