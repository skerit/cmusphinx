#include <stdio.h>
#include <string.h>
#include <time.h>

#include "pocketsphinx_internal.h"
#include "bptbl.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	ps_decoder_t *ps;
	cmd_ln_t *config;
	bptbl_t *bptbl;

	/* Get the API to initialize a bunch of stuff for us (but not the search). */
	config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
			     "-dict", TESTDATADIR "/hub4.5000.dic",
			     "-fwdtree", "no",
			     "-fwdflat", "no",
			     "-bestpath", "no", NULL);
	ps = ps_init(config);
	bptbl = bptbl_init(ps->d2p, 10, 10);

	bptbl_free(bptbl);
	ps_free(ps);
}
