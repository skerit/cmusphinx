/**
 * @file test_decoder_factory.c
 * @author dhuggins
 */

#include <multisphinx/search_factory.h>

#include "test_macros.h"

int main(int argc, char **argv)
{
    search_factory_t *dcf;
    search_t *fwdtree;

    dcf = search_factory_init(
            "-lm", TESTDATADIR "/bn10000.3g.arpa",
            "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
            "-dict", TESTDATADIR "/bn10000.dic",
            "-samprate", "11025", NULL);
    TEST_ASSERT(dcf != NULL);

    fwdtree = search_factory_create(dcf, "fwdtree", NULL);
    TEST_ASSERT(fwdtree != NULL);

    return 0;
}
