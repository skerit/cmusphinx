/**
 * @file test_state_align.c
 * @author dhuggins
 */

#include <multisphinx/search_factory.h>
#include <multisphinx/featbuf.h>
#include <multisphinx/state_align_search.h>

#include "test_macros.h"

int main(int argc, char *argv[])
{
    search_factory_t *dcf;
    alignment_t *al;
    featbuf_t *fb;
    search_t *sas;
    FILE *infh;
    size_t nread;
    int16 buf[2048];

    dcf = search_factory_init("-lm", TESTDATADIR "/bn10000.3g.arpa", "-hmm",
            TESTDATADIR "/hub4wsj_sc_8k", "-dict", TESTDATADIR "/bn10000.dic",
            "-samprate", "11025", NULL);
    TEST_ASSERT(dcf != NULL);
    sas = search_factory_create(dcf, "state_align", NULL);
    TEST_ASSERT(sas != NULL);

    al = alignment_init(search_factory_d2p(dcf));
    alignment_add_words(al, "<s>", "A", "LONG", "DISCUSSION", "ABOUT", "HOW", "WE",
            "WANT", "TO", "MAKE", "IT", "EASIER", "FOR", "PEOPLE", "TO",
            "BLEEP", "THINGS", "OUT", "</s>", NULL);
    alignment_populate_ci(al);
    state_align_search_set_alignment(sas, al);
    search_run(sas);

    fb = search_factory_featbuf(dcf);
    TEST_ASSERT(fb != NULL);
    infh = fopen(TESTDATADIR "/chan3.raw", "rb");
    TEST_ASSERT(infh != NULL);
    featbuf_producer_start_utt(fb, NULL);
    nread = 0;
    while (nread < 77175)
    {
        size_t n = fread(buf, 2, 2048, infh);
        if (n == 0)
            break;
        featbuf_producer_process_raw(fb, buf, n, FALSE);
        nread += n;
    }
    fclose(infh);
    featbuf_producer_end_utt(fb);
    featbuf_producer_shutdown(fb);

    {
        int32 score;
        char const *hyp = search_hyp(sas, &score);
        printf("%s (%d)\n", hyp, score);
    }

    {
        int32 score;
        seg_iter_t *itor = search_seg_iter(sas, &score);
        for (; itor; itor = seg_iter_next(itor)) {
            char const *word = seg_iter_word(itor);
            int sf, ef;
            seg_iter_times(itor, &sf, &ef);
            printf("%s %d:%d\n", word, sf, ef);
        }
    }
    search_free(sas);
    search_factory_free(dcf);

    return 0;
}
