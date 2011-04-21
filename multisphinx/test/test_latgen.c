#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/feat.h>

#include <multisphinx/pocketsphinx.h>
#include <multisphinx/pocketsphinx_internal.h>
#include <multisphinx/fwdtree_search.h>
#include <multisphinx/fwdflat_search.h>
#include <multisphinx/latgen_search.h>

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
	featbuf_t *fb;
	search_t *fwdtree, *fwdflat, *latgen;
	FILE *mfcfh;
	size_t nsamp;
	char const *hyp;
	int32 score;
	mfcc_t buf[13];

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/bn10000.3g.homos.arpa",
			     "-dict", TESTDATADIR "/bn10000.homos.dic",
			     "-arcdumpdir", ".",
			     "-maxwpf", "50",
			     "-latsize", "512",
			     NULL);
	ps_init_defaults(config);
	fb = featbuf_init(config);

	/* Create acoustic model and search. */
	lmath = logmath_init(cmd_ln_float32_r(config, "-logbase"),
			     0, FALSE);
	acmod = acmod_init(config, lmath, fb);
	acmod2 = acmod_copy(acmod);
	mdef = bin_mdef_read(config, cmd_ln_str_r(config, "-mdef"));
	dict = dict_init(config, mdef);
	d2p = dict2pid_build(mdef, dict);
	fwdtree = fwdtree_search_init(config, acmod, dict, d2p);
	fwdflat = fwdflat_search_init(config, acmod2, dict, d2p,
				      ps_search_lmset(fwdtree));
	latgen = latgen_init(config, d2p, ps_search_lmset(fwdtree));
	ps_search_link(fwdtree, fwdflat, "fwdtree", FALSE);
	ps_search_link(fwdflat, latgen, "fwdflat", TRUE);

	/* Launch search threads. */
	ps_search_run(fwdtree);
	ps_search_run(fwdflat);
	ps_search_run(latgen);

	/* Feed them a bunch of data. */
	if ((mfcfh = fopen(TESTDATADIR "/050c0103.mfc", "rb")) == NULL) {
		E_FATAL_SYSTEM("Failed to open "TESTDATADIR"/050c0103.mfc");
		return 1;
	}
	fread(buf, 4, 1, mfcfh);
	featbuf_producer_start_utt(fb, "050c0103");
	while ((nsamp = fread(buf, 4, 13, mfcfh)) > 0)  {
		mfcc_t *bptr = buf;
		featbuf_producer_process_cep(fb, &bptr, 1, FALSE);
	}
	fclose(mfcfh);

	/* This will wait for search to complete. */
	E_INFO("Waiting for end of utt\n");
	featbuf_producer_end_utt(fb);
	E_INFO("Done waiting\n");

	/* Retrieve the hypothesis from the search thread. */
	hyp = ps_search_hyp(fwdflat, &score);
	E_INFO("hyp: %s (%d)\n", hyp, score);

	/* Reap the search thread. */
	E_INFO("Reaping the search threads\n");
	featbuf_producer_shutdown(fb);
	ps_search_wait(fwdtree);
	ps_search_wait(fwdflat);
	ps_search_wait(latgen);
	E_INFO("Done reaping\n");
	ps_search_free(fwdtree);
	ps_search_free(fwdflat);
	ps_search_free(latgen);
	acmod_free(acmod);
	acmod_free(acmod2);
	featbuf_free(fb);

	/* Clean everything else up. */
	dict_free(dict);
	bin_mdef_free(mdef);
	dict2pid_free(d2p);
	logmath_free(lmath);
	cmd_ln_free_r(config);

	return 0;
}
