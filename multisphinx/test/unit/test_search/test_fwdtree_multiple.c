#include <pocketsphinx.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/feat.h>

#include "pocketsphinx_internal.h"
#include "fwdtree_search.h"
#include "fwdflat_search.h"
#include "test_macros.h"

static ps_search_t *
init_search(cmd_ln_t *config, acmod_t *acmod,
	    dict_t *dict, dict2pid_t *d2p,
	    char const *lmfile)
{
	acmod_t *acmod2;
	ps_search_t *search;

	acmod2 = acmod_copy(acmod);
	cmd_ln_set_str_r(config, "-lm", lmfile);
	search = fwdtree_search_init(config, acmod2, dict, d2p);
	acmod_free(acmod2);

	return search;
}

int
main(int argc, char *argv[])
{
	acmod_t *acmod;
	logmath_t *lmath;
	cmd_ln_t *config;
	featbuf_t *fb;
	ps_search_t *searches[4];
	mfcc_t ***feat;
	int nfr, i;
	char const *hyp;
	int32 score;
	dict2pid_t *d2p;
	dict_t *dict;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-dict", TESTDATADIR "/cmu07a.dic",
			     NULL);
	ps_init_defaults(config);
	fb = featbuf_init(config);
	lmath = logmath_init(cmd_ln_float32_r(config, "-logbase"),
			     0, FALSE);
	acmod = acmod_init(config, lmath, fb);

	/* Construct searches for 5, 10, 20, and 40k vocab. */
	dict = dict_init(config, acmod->mdef);
	d2p = dict2pid_build(acmod->mdef, dict);
	searches[0] = init_search(config, acmod, dict, d2p,
				  TESTDATADIR "/bn.5000.arpa.DMP");
	searches[1] = init_search(config, acmod, dict, d2p,
				  TESTDATADIR "/bn.10000.arpa.DMP");
	searches[2] = init_search(config, acmod, dict, d2p,
				  TESTDATADIR "/bn.20000.arpa.DMP");
	searches[3] = init_search(config, acmod, dict, d2p,
				  TESTDATADIR "/bn.40000.arpa.DMP");
	
	/* Launch search threads. */
	for (i = 0; i < 4; ++i)
		ps_search_run(searches[i]);

	/* Feed them a bunch of data. */
	nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
			      ".mfc", 0, -1, NULL, -1);
	feat = feat_array_alloc(acmod->fcb, nfr);
	if ((nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
				   ".mfc", 0, -1, feat, -1)) < 0)
		E_FATAL("Failed to read mfc file\n");
	acmod_free(acmod);

	featbuf_start_utt(fb);
	for (i = 0; i < nfr; ++i)
		featbuf_process_feat(fb, feat[i]);

	/* Wait for searches to complete. */
	E_INFO("Waiting for end of utt\n");
	featbuf_end_utt(fb, -1);
	E_INFO("Done waiting\n");

	/* Retrieve hypotheses and timing information. */
	for (i = 0; i < 4; ++i) {
		hyp = ps_search_hyp(searches[i], &score);
		E_INFO("hyp %d: %s (%d)\n", i, hyp, score);
	}

	/* Reap the search threads. */
	E_INFO("Reaping the search threads\n");
	featbuf_shutdown(fb);
	E_INFO("Done reaping\n");
	for (i = 0; i < 4; ++i)
		ps_search_free(searches[i]);
	featbuf_free(fb);

	/* Clean everything else up. */
	feat_array_free(feat);
	logmath_free(lmath);
	cmd_ln_free_r(config);
	dict_free(dict);
	dict2pid_free(d2p);

	return 0;
}
