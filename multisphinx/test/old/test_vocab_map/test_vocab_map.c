/**
 * @file test_vocab_map.c
 * @brief Test the vocabulary mapping object.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include <multisphinx/pocketsphinx_internal.h>
#include <multisphinx/vocab_map.h>
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	vocab_map_t *vm;
	cmd_ln_t *config;
	bin_mdef_t *mdef;
	FILE *vmfh;
	dict_t *dict;
	int32 const *wids;
	int32 nmapped;
	int32 dhehr;

	config = cmd_ln_init(NULL, ps_args(), FALSE,
			     "-dict", TESTDATADIR "/bn10000.homos.dic",
			     NULL);
	mdef = bin_mdef_read(config, TESTDATADIR "/hub4wsj_sc_8k/mdef");
	dict = dict_init(config, mdef);
	TEST_ASSERT(vm = vocab_map_init(dict));
	vmfh = fopen(TESTDATADIR "/bn10000.homos", "r");
	TEST_ASSERT(0 == vocab_map_read(vm, vmfh));
	fclose(vmfh);

	vmfh = fopen("tmp.bn10000.homos", "w");
	vocab_map_write(vm, vmfh);
	fclose(vmfh);
	TEST_ASSERT(0 == system("LANG=C sort '"TESTDATADIR"/bn10000.homos' "
				"| diff -w -u - tmp.bn10000.homos"));

	/* FIXME: Actually we should probably test every single word. */
	dhehr = dict_wordid(dict, "_DH_EH_R");
	wids = vocab_map_unmap(vm, dhehr, &nmapped);
	TEST_ASSERT(wids);
	TEST_ASSERT(nmapped == 4);
	TEST_ASSERT(wids[0] == dict_wordid(dict, "THERE."));
	TEST_ASSERT(dhehr == vocab_map_map(vm, wids[0]));
	TEST_ASSERT(wids[1] == dict_wordid(dict, "THERE"));
	TEST_ASSERT(dhehr == vocab_map_map(vm, wids[1]));
	TEST_ASSERT(wids[2] == dict_wordid(dict, "THEY'RE"));
	TEST_ASSERT(dhehr == vocab_map_map(vm, wids[2]));
	TEST_ASSERT(wids[3] == dict_wordid(dict, "THEIR"));
	TEST_ASSERT(dhehr == vocab_map_map(vm, wids[3]));

	unlink("tmp.bn10000.homos");
	dict_free(dict);
	vocab_map_free(vm);
	bin_mdef_free(mdef);
	cmd_ln_free_r(config);

	return 0;
}
