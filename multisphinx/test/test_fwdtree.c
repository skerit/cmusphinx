#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/feat.h>

#include <multisphinx/pocketsphinx.h>
#include <multisphinx/fwdtree_search.h>
#include <multisphinx/acmod.h>
#include <multisphinx/ps_search.h>
#include <multisphinx/featbuf.h>

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
	featbuf_t *fb;
	ps_search_t *fwdtree;
	mfcc_t ***feat;
	int nfr, i;
	char const *hyp;
	int32 score;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/bn10000.3g.arpa",
			     "-dict", TESTDATADIR "/bn10000.dic",
			     "-samprate", "11025",
			     NULL);
	ps_init_defaults(config);
	fb = featbuf_init(config);

	/* Create acoustic model and search. */
	lmath = logmath_init(cmd_ln_float32_r(config, "-logbase"),
			     0, FALSE);
	acmod = acmod_init(config, lmath, fb);
	mdef = bin_mdef_read(config, cmd_ln_str_r(config, "-mdef"));
	dict = dict_init(config, mdef);
	d2p = dict2pid_build(mdef, dict);
	fwdtree = fwdtree_search_init(config, acmod, dict, d2p);

	/* Launch a search thread. */
	ps_search_run(fwdtree);

	/* Feed it a bunch of data. */
	nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
			      ".mfc", 0, -1, NULL, -1);
	feat = feat_array_alloc(acmod->fcb, nfr);
	if ((nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
				   ".mfc", 0, -1, feat, -1)) < 0)
		E_FATAL("Failed to read mfc file\n");
	featbuf_producer_start_utt(fb, "chan3");
	for (i = 0; i < nfr; ++i)
		featbuf_producer_process_feat(fb, feat[i]);

	/* This will wait for search to complete. */
	printf("Waiting for end of utt\n");
	featbuf_producer_end_utt(fb);
	printf("Done waiting\n");

	/* Retrieve the hypothesis from the search thread. */
	hyp = ps_search_hyp(fwdtree, &score);
	printf("hyp: %s (%d)\n", hyp, score);

	/* Reap the search thread. */
	printf("Reaping the search thread\n");
	featbuf_producer_shutdown(fb);
	ps_search_wait(fwdtree);
	printf("Done reaping\n");
	ps_search_free(fwdtree);
	acmod_free(acmod);
	featbuf_free(fb);

	/* Clean everything else up. */
	dict_free(dict);
	bin_mdef_free(mdef);
	dict2pid_free(d2p);
	feat_array_free(feat);
	logmath_free(lmath);
	cmd_ln_free_r(config);

	return 0;
}
