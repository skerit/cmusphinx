/**
 * @file test_decoder_factory.c
 * @author dhuggins
 */

#include <sphinxbase/feat.h>

#include <multisphinx/search_factory.h>
#include <multisphinx/acmod.h>

#include "test_macros.h"

int main(int argc, char **argv)
{
    search_factory_t *dcf;
    featbuf_t *fb;
    search_t *fwdtree;
    acmod_t *acmod;
    int nfr, i;
    mfcc_t ***feat;

    dcf = search_factory_init(
            "-lm", TESTDATADIR "/bn10000.3g.arpa",
            "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
            "-dict", TESTDATADIR "/bn10000.dic",
            "-samprate", "11025", NULL);
    TEST_ASSERT(dcf != NULL);

    fwdtree = search_factory_create(dcf, "fwdtree", NULL);
    TEST_ASSERT(fwdtree != NULL);

    search_run(fwdtree);

    fb = search_factory_featbuf(dcf);
    TEST_ASSERT(fb!= NULL);
    acmod = search_factory_acmod(dcf);

    /* Feed it a bunch of data. */
    nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
                          ".mfc", 0, -1, NULL, -1);
    feat = feat_array_alloc(acmod->fcb, nfr);
    if ((nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR,
                               ".mfc", 0, -1, feat, -1)) < 0) {
            E_ERROR("Failed to read mfc file\n");
            return 1;
    }

    featbuf_producer_start_utt(fb, NULL);
    for (i = 0; i < nfr; ++i)
            featbuf_producer_process_feat(fb, feat[i]);

    /* This will wait for search to complete. */
    printf("Waiting for end of utt\n");
    featbuf_producer_end_utt(fb);
    printf("Done waiting\n");

    /* Retrieve the hypothesis from the search thread. */
    //hyp = ps_search_hyp(fwdflat, &score);
    //printf("hyp: %s (%d)\n", hyp, score);

    /* Reap the search thread. */
    E_INFO("Reaping the search thread\n");
    featbuf_producer_shutdown(fb);
    printf("Done reaping\n");
    search_free(fwdtree);
    acmod_free(acmod);
    featbuf_free(fb);


    return 0;
}
