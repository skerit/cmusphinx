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
    char *hyp1, *hyp2;
    int32 score;

    dcf = search_factory_init("-lm", TESTDATADIR "/bn10000.3g.arpa", "-hmm",
            TESTDATADIR "/hub4wsj_sc_8k", "-dict", TESTDATADIR "/bn10000.dic",
            "-samprate", "11025", NULL);
    TEST_ASSERT(dcf != NULL);

    fwdtree = search_factory_create(dcf, "fwdtree", NULL);
    TEST_ASSERT(fwdtree != NULL);

    search_run(fwdtree);

    fb = search_factory_featbuf(dcf);
    TEST_ASSERT(fb!= NULL);
    acmod = search_factory_acmod(dcf);

    /* Feed it a bunch of data. */
    nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR, ".mfc", 0, -1,
            NULL, -1);
    feat = feat_array_alloc(acmod->fcb, nfr);
    if ((nfr = feat_s2mfc2feat(acmod->fcb, "chan3", TESTDATADIR, ".mfc", 0, -1,
            feat, -1)) < 0)
    {
        E_ERROR("Failed to read mfc file\n");
        return 1;
    }

    featbuf_producer_start_utt(fb, NULL);
    for (i = 0; i < 500; ++i)
        featbuf_producer_process_feat(fb, feat[i]);

    /* This will wait for search to complete. */
    featbuf_producer_end_utt(fb);
    /* Retrieve the hypothesis from the search thread. */
    hyp1 = ckd_salloc(search_hyp(fwdtree, &score));
    E_INFO("hyp: %s (%d)\n", hyp1, score);

    /* Reap the search thread. */
    featbuf_producer_shutdown(fb);
    search_free(fwdtree);

    /* Now verify that overriding (and also creating new searches) works */
    fwdtree = search_factory_create(dcf, "fwdtree", "-lm", TESTDATADIR "/bn10000.3g.homos.arpa", NULL);
    search_run(fwdtree);
    featbuf_producer_start_utt(fb, NULL);
    for (i = 0; i < 500; ++i)
        featbuf_producer_process_feat(fb, feat[i]);
    featbuf_producer_end_utt(fb);

    /* Retrieve the hypothesis from the search thread. */
    hyp2 = ckd_salloc(search_hyp(fwdtree, &score));
    E_INFO("hyp 2: %s (%d)\n", hyp2, score);
    TEST_ASSERT(0 != strcmp(hyp1, hyp2));
    featbuf_producer_shutdown(fb);
    search_free(fwdtree);
    search_factory_free(dcf);
    ckd_free(hyp1);
    ckd_free(hyp2);

    return 0;
}
