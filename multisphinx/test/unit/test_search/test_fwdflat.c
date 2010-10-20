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
	ps_decoder_t *ps;
	cmd_ln_t *config;
	acmod_t *acmod;
	ps_search_t *fwdflat;
	mfcc_t ***feat;
	int nfr, i;
	char const *hyp;
	int32 score;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     "-fwdtree", "no",
			     "-fwdflat", "no",
			     "-bestpath", "no", NULL);

	/* Get the API to initialize a bunch of stuff for us (but not the search). */
	ps = ps_init(config);
	acmod = ps->acmod;
	cmd_ln_set_str_r(config, "-lm", TESTDATADIR "/hub4.5000.DMP");
	fwdflat = fwdflat_search_init(config, acmod, ps->dict, ps->d2p, NULL);

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
	ps_free(ps);

	return 0;
}
