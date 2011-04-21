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
	acmod_t *acmod2, *acmod3;
	ps_search_t *search, *search2;

	acmod2 = acmod_copy(acmod);
	cmd_ln_set_str_r(config, "-lm", lmfile);
	search = fwdtree_search_init(config, acmod2, dict, d2p);
	acmod3 = acmod_copy(acmod);
	search2 = fwdflat_search_init(config, acmod3, dict, d2p,
				      ps_search_output_arcs(search),
				      fwdtree_search_lmset(search));
	acmod_free(acmod2);
	acmod_free(acmod3);
	
	/* Launch search threads. */
	ps_search_run(search);
	ps_search_run(search2);

	return search2;
}

int
main(int argc, char *argv[])
{
	acmod_t *acmod;
	logmath_t *lmath;
	cmd_ln_t *config;
	featbuf_t *fb;
	ps_search_t *searches[4];
	char const *hyp;
	int32 score;
	dict2pid_t *d2p;
	dict_t *dict;
	FILE *rawfh;
	FILE *hypfh;
	int16 buf[2048];
	size_t nsamp, tsamp;
	ptmr_t total;
	int i;

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
				  TESTDATADIR "/bn.2000.arpa.DMP");
	searches[1] = init_search(config, acmod, dict, d2p,
				  TESTDATADIR "/bn.5000.arpa.DMP");
	searches[2] = init_search(config, acmod, dict, d2p,
				  TESTDATADIR "/bn.10000.arpa.DMP");
	searches[3] = init_search(config, acmod, dict, d2p,
				  TESTDATADIR "/bn.20000.arpa.DMP");
	acmod_free(acmod);

	/* Feed them a bunch of data. */
	if ((rawfh = fopen(TESTDATADIR "/i960711p.raw", "rb")) == NULL) {
		E_FATAL_SYSTEM("Failed to open "TESTDATADIR"/i960711p.raw");
		return 1;
	}
	ptmr_init(&total);
	ptmr_start(&total);
	featbuf_producer_start_utt(fb, "i960711p");
	tsamp = 0;
	while ((nsamp = fread(buf, 2, 2048, rawfh)) > 0) {
		tsamp += nsamp;
		featbuf_producer_process_raw(fb, buf, nsamp, FALSE);
	}
	fclose(rawfh);
	/* Wait for searches to complete. */
	E_INFO("Waiting for end of utt\n");
	featbuf_producer_end_utt(fb, -1);
	E_INFO("Done waiting\n");
	ptmr_stop(&total);
	E_INFO("%f secs of speech, %f elapsed, %f CPU\n",
	       (double)tsamp / 16000,
	       total.t_elapsed, total.t_cpu);

	/* Retrieve hypotheses and timing information. */
	if ((hypfh = fopen("test_fwdtree_multiple.hyp", "w")) == NULL) {
		E_FATAL_SYSTEM("Failed to open test_fwdtree_multiple.hyp");
		return 1;
	}
	for (i = 0; i < 4; ++i) {
		hyp = ps_search_hyp(searches[i], &score);
		E_INFO("hyp %d: %s (%d)\n", i, hyp, score);
		fprintf(hypfh, "%s (i960711p_%d %d)\n", hyp, i, score);
	}
	fclose(hypfh);

	/* Reap the search threads. */
	E_INFO("Reaping the search threads\n");
	featbuf_producer_shutdown(fb);
	for (i = 0; i < 4; ++i) {
		ps_search_wait(searches[i]);
		ps_search_free(searches[i]);
	}
	E_INFO("Done reaping\n");
	featbuf_free(fb);

	/* Clean everything else up. */
	logmath_free(lmath);
	cmd_ln_free_r(config);
	dict_free(dict);
	dict2pid_free(d2p);

	return 0;
}
