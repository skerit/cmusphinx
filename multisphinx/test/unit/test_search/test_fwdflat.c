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
	fwdtree_search_t *fwdtree;
	mfcc_t **feat;
	int nfr;

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
	fwdtree = fwdtree_search_init(config, acmod, ps->dict, ps->dict2pid);

	if ((nfr = feat_s2mfc2feat(acmod->feat, "chan3", TESTDATADIR,
				   ".mfc", 0, -1, &feat, -1) < 0))
		E_FATAL("Failed to read MFC file " TESTDATADIR "/chan3.mfc\n");
	ps_search_start(fwdtree);
	for (i = 0; i < nfr; ++i) {
		acmod_process_feat(acmod, feat + i);
		while (acmod->n_feat_frame > 0) {
			ps_search_step(fwdtree, acmod->output_frame);
			acmod_advance(acmod);
		}
	}
	ps_search_finish(fwdtree);

	ps_search_free(fwdtree);
	ps_free(ps);

	return 0;
}
