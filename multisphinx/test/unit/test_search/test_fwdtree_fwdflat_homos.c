#include <pocketsphinx.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sphinxbase/feat.h>

#include <multisphinx/fwdtree_search.h>
#include <multisphinx/fwdflat_search.h>

#include "pocketsphinx_internal.h"
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
	ps_search_t *fwdtree, *fwdflat;
	vocab_map_t *vm;
	FILE *rawfh;
	int16 buf[2048];
	size_t nsamp;
	char const *hyp;
	int32 score;

	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-lm", TESTDATADIR "/bn10000.3g.homos.arpa",
			     "-dict", TESTDATADIR "/bn10000.homos.dic",
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
	cmd_ln_set_str_r(config, "-lm", TESTDATADIR "/bn10000.3g.arpa");
	fwdflat = fwdflat_search_init(config, acmod2, dict, d2p,
				      ps_search_output_arcs(fwdtree),
				      NULL);
	vm = vocab_map_init(dict);
	rawfh = fopen(TESTDATADIR "/bn10000.homos", "r");
	vocab_map_read(vm, rawfh);
	fclose(rawfh);
	fwdflat_search_set_vocab_map(fwdflat, vm);

	/* Launch search threads. */
	ps_search_run(fwdtree);
	ps_search_run(fwdflat);

	/* Feed them a bunch of data. */
	if ((rawfh = fopen(TESTDATADIR "/i960711p.raw", "rb")) == NULL) {
		E_FATAL_SYSTEM("Failed to open "TESTDATADIR"/i960711p.raw");
		return 1;
	}
	featbuf_producer_start_utt(fb, "i960711p");
	while ((nsamp = fread(buf, 2, 2048, rawfh)) > 0)
		featbuf_producer_process_raw(fb, buf, nsamp, FALSE);
	fclose(rawfh);

	/* This will wait for search to complete. */
	E_INFO("Waiting for end of utt\n");
	featbuf_producer_end_utt(fb, -1);
	E_INFO("Done waiting\n");

	/* Retrieve the hypothesis from the search thread. */
	hyp = ps_search_hyp(fwdflat, &score);
	E_INFO("hyp: %s (%d)\n", hyp, score);

	/* Reap the search thread. */
	E_INFO("Reaping the search threads\n");
	featbuf_producer_shutdown(fb);
	ps_search_wait(fwdtree);
	ps_search_wait(fwdflat);
	E_INFO("Done reaping\n");
	ps_search_free(fwdflat);
	acmod_free(acmod);
	featbuf_free(fb);

	/* Clean everything else up. */
	dict_free(dict);
	bin_mdef_free(mdef);
	dict2pid_free(d2p);
	logmath_free(lmath);
	cmd_ln_free_r(config);

	return 0;
}
