/**
 * @file test_decoder_factory.c
 * @author dhuggins
 */

#include <multisphinx/decoder_factory.h>

#include "test_macros.h"

int main(int argc, char **argv)
{
    search_factory_t *dcf;
    config_t *config;

    config = search_factory_config(
            "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
            "-lm", TESTDATADIR "/bn10000.3g.arpa",
            "-dict", TESTDATADIR "/bn10000.dic",
            "-samprate", "11025", NULL);
    TEST_ASSERT(config != NULL);
    dcf = search_factory_init(config);
    TEST_ASSERT(dcf != NULL);

    return 0;
}
