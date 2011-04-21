/**
 * @file test_partial_results.c
 * @author dhuggins
 */

#include <sphinxbase/feat.h>

#include <multisphinx/search_factory.h>

#include <sys/time.h>
#include <unistd.h>

#include "test_macros.h"

static struct timeval utt_start;

static int
search_cb(search_t *search, search_event_t *evt, void *udata)
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    E_INFO("%s: time delta %f ",
            search_name(search),
           ((double)tv.tv_sec + (double)tv.tv_usec / 1000000)
           - ((double)utt_start.tv_sec + (double)utt_start.tv_usec / 1000000));
    switch (evt->event) {
    case SEARCH_PARTIAL_RESULT: {
        int32 score;
        char const *hyp = search_hyp(search, &score);
        E_INFOCONT("partial: %s (%d)\n", hyp, score);
        break;
    }
    case SEARCH_START_UTT:
        E_INFOCONT("start\n");
        break;
    case SEARCH_END_UTT:
        E_INFOCONT("end\n");
        break;
    case SEARCH_FINAL_RESULT: {
        int32 score;
        char const *hyp = search_hyp(search, &score);
        E_INFOCONT("full: %s (%d)\n", hyp, score);
        break;
    }
    }
    return 0;
}

int main(int argc, char **argv)
{
    search_factory_t *dcf;
    featbuf_t *fb;
    search_t *fwdtree, *fwdflat;
    FILE *infh;
    size_t nread;
    int16 buf[2048];

    dcf = search_factory_init("-lm", TESTDATADIR "/bn10000.3g.arpa", "-hmm",
            TESTDATADIR "/hub4wsj_sc_8k", "-dict", TESTDATADIR "/bn10000.dic",
            "-samprate", "11025", NULL);
    TEST_ASSERT(dcf != NULL);

    fwdtree = search_factory_create(dcf, "fwdtree", NULL);
    fwdflat = search_factory_create(dcf, "fwdflat", NULL);
    search_link(fwdtree, fwdflat, "fwdtree", FALSE);

    search_set_cb(fwdtree, search_cb, NULL);
    search_set_cb(fwdflat, search_cb, NULL);
    search_run(fwdtree);
    search_run(fwdflat);

    fb = search_factory_featbuf(dcf);
    TEST_ASSERT(fb != NULL);

    /* Feed it a bunch of data. */
    gettimeofday(&utt_start, NULL);
    infh = fopen(TESTDATADIR "/chan3.raw", "rb");
    TEST_ASSERT(infh != NULL);
    featbuf_producer_start_utt(fb, NULL);
    nread = 0;
    /* Magic numbers!  Process 7 seconds of audio, and sleep for 180ms (~= 2048 / 11025) */
    while (nread < 77175) {
        size_t n = fread(buf, 2, 2048, infh);
        if (n == 0)
            break;
        featbuf_producer_process_raw(fb, buf, n, FALSE);
        usleep(180000);
        nread += n;
    }
    fclose(infh);
    featbuf_producer_end_utt(fb);
    featbuf_producer_shutdown(fb);
    search_free(fwdtree);
    search_free(fwdflat);
    search_factory_free(dcf);

    return 0;
}
