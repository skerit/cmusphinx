#include <multisphinx/arc_buffer.h>
#include <multisphinx/fwdflat_search.h>

#include "test_macros.h"

static logmath_t *lmath;
static bin_mdef_t *mdef;
static dict2pid_t *d2p;
static dict_t *dict;
static cmd_ln_t *config;
static ngram_model_t *lm;

static void test_arcbuf(arc_buffer_t *arcs)
{
    int fi, i, next_sf;
    bptbl_t *bptbl;
    arc_t *a, *aa;
    bpidx_t bp;
    bp_t bpe;

    bptbl = arc_buffer_input_bptbl(arcs);

    /* Enter a bunch of initial bps (like silence) */
    fi = bptbl_push_frame(bptbl, NO_BP);
    bp = bptbl_enter(bptbl, 42, NO_BP, 1, 0);
    fi = bptbl_push_frame(bptbl, NO_BP);
    bp = bptbl_enter(bptbl, 42, NO_BP, 2, 0);
    fi = bptbl_push_frame(bptbl, NO_BP);
    bp = bptbl_enter(bptbl, 42, NO_BP, 3, 0);

    /* Enter a couple of words pointing back to the silences. */
    fi = bptbl_push_frame(bptbl, NO_BP);
    bp = bptbl_enter(bptbl, 69, 1, 4, 0);
    bp = bptbl_enter(bptbl, 69, 1, 5, 0);

    /* Garbage collect some things. */
    fi = bptbl_push_frame(bptbl, 2);

    /* Add some more words. */
    fi = bptbl_push_frame(bptbl, 2);
    bp = bptbl_enter(bptbl, 999, 4, 5, 0);

    /* Garbage collect some things. */
    fi = bptbl_push_frame(bptbl, 5);
    bptbl_dump(bptbl);
    next_sf = bptbl_active_sf(bptbl);
    E_INFO("next_sf %d\n", next_sf);
    arc_buffer_producer_sweep(arcs, FALSE);

    /* Now add a bunch of stuff to see what happens. */
    for (i = 0; i < 6; ++i)
    {
        bp = bptbl_enter(bptbl, 42, 5, 6 + i, 0);
    }
    fi = bptbl_push_frame(bptbl, 9);
    for (i = 0; i < 3; ++i)
    {
        bp = bptbl_enter(bptbl, 69, 6, 12 + i, 0);
    }
    fi = bptbl_push_frame(bptbl, 12);
    bptbl_dump(bptbl);
    bptbl_get_bp(bptbl, bptbl->oldest_bp, &bpe);
    next_sf = bpe.frame + 1;
    E_INFO("next_sf %d\n", next_sf);
    arc_buffer_producer_sweep(arcs, FALSE);

    for (i = 0; i < 3; ++i)
    {
        bp = bptbl_enter(bptbl, 420, 6, 39 + i, 0);
    }
    bptbl_finalize(bptbl);
    bptbl_dump(bptbl);
    bptbl_get_bp(bptbl, bptbl->oldest_bp, &bpe);
    next_sf = bpe.frame + 1;
    E_INFO("next_sf %d\n", next_sf);
    arc_buffer_producer_sweep(arcs, FALSE);
    arc_buffer_dump(arcs, dict);

    aa = a = arc_buffer_iter(arcs, 2);
    TEST_ASSERT(a->wid == 69);
    a = arc_buffer_iter_next(arcs, a);
    a = arc_buffer_iter_next(arcs, a);
    E_INFO("a was %p, is now %p (%d bytes) arcs[4] = %p\n",
            aa, a, (char *)a - (char *)aa,
            arc_buffer_iter(arcs, 4));
    /* There should be two arcs exiting frame 2 and none in frame 3. */
    TEST_ASSERT(a == arc_buffer_iter(arcs, 4));

    for (a = arc_buffer_iter(arcs, 6);
            a != arc_buffer_iter(arcs, 8);
            a = arc_buffer_iter_next(arcs, a))
    {
        TEST_ASSERT(a->wid == 42 || a->wid == 420);
        TEST_ASSERT(a->src >= 6 && a->src < 8);
    }

    arc_buffer_dump(arcs, dict);
    arc_buffer_release(arcs, 6);
    for (a = arc_buffer_iter(arcs, 6);
            a != arc_buffer_iter(arcs, 8);
            a = arc_buffer_iter_next(arcs, a))
    {
        printf("%d %d %d\n", a->wid, a->src, a->dest);
        TEST_ASSERT(a->wid == 42 || a->wid == 420);
        TEST_ASSERT(a->src >= 6 && a->src < 8);
    }
    arc_buffer_dump(arcs, dict);
}

int main(int argc, char *argv[])
{
    decoder_factory_t *dcf;
    arc_buffer_t *arcs;
    bptbl_t *bptbl;

    /* Get the API to initialize a bunch of stuff for us (but not the search). */
    config = decoder_factory_config(
            "-hmm", TESTDATADIR "/hub4wsj_sc_8k",
            "-lm", TESTDATADIR "/hub4.5000.DMP",
            "-dict", TESTDATADIR "/hub4.5000.dic",
            NULL);
    dcf = decoder_factory_init(config);
    mdef = acmod_mdef(decoder_factory_acmod(dcf));
    d2p = decoder_factory_d2p(dcf);
    lm = decoder_factory_lm(dcf);

    bptbl = bptbl_init("test", d2p, 10, 10);
    arcs = arc_buffer_init("noscore", bptbl, lm, FALSE);
    test_arcbuf(arcs);
    arc_buffer_free(arcs);
    bptbl_free(bptbl);

    bptbl = bptbl_init("test", d2p, 10, 10);
    arcs = arc_buffer_init("score", bptbl, lm, TRUE);
    test_arcbuf(arcs);
    arc_buffer_free(arcs);
    bptbl_free(bptbl);

    decoder_factory_free(dcf);

    return 0;
}
