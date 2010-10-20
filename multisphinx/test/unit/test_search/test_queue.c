#include <pocketsphinx.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/feat.h>

#include "pocketsphinx_internal.h"
#include "queue_search.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	ps_decoder_t *ps;
	cmd_ln_t *config;
	acmod_t *acmod;
	fwdtree_search_t *fwdtree;
	fwdflat_search_t *fwdflat;
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

	/* Initialize search modules. */
	cmd_ln_set_str_r(config, "-lm", TESTDATADIR "/hub4.5000.DMP");
	fwdtree = fwdtree_search_init(config, acmod, ps->dict, ps->dict2pid);
	fwdflat = fwdflat_search_init(config, acmod, ps->dict, ps->dict2pid);
	searchq = queue_search_init(config);

	if ((nfr = feat_s2mfc2feat(acmod->feat, "chan3", TESTDATADIR,
				   ".mfc", 0, -1, &feat, -1) < 0))
		E_FATAL("Failed to read MFC file " TESTDATADIR "/chan3.mfc\n");
	ps_search_start(searchq);
	for (i = 0; i < nfr; ++i) {
		acmod_process_feat(acmod, feat + i);
		while (acmod->n_feat_frame > 0) {
			ps_search_step(searchq, acmod->output_frame);
			acmod_advance(acmod);
		}
	}
	ps_search_finish(searchq);
	ps_search_free(searchq);
	ps_free(ps);

	return 0;
}
