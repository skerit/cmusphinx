/**
 * @file test_vocab_map.c
 * @brief Test the vocabulary mapping object.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include <pocketsphinx.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "pocketsphinx_internal.h"
#include "vocab_map.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	vocab_map_t *vm;
	cmd_ln_t *config;
	bin_mdef_t *mdef;
	FILE *vmfh;
	dict_t *dict;

	config = cmd_ln_init(NULL, ps_args(), FALSE,
			     "-dict", TESTDATADIR "/bn10000.homos.dic",
			     NULL);
	mdef = bin_mdef_read(config, TESTDATADIR "/hub4wsj_sc_8k/mdef");
	dict = dict_init(config, mdef);
	TEST_ASSERT(vm = vocab_map_init(dict));
	vmfh = fopen(TESTDATADIR "/bn10000.homos", "r");
	TEST_ASSERT(0 == vocab_map_read(vm, vmfh));
	fclose(vmfh);

	dict_free(dict);
	vocab_map_free(vm);
	bin_mdef_free(mdef);
	cmd_ln_free_r(config);

	return 0;
}
