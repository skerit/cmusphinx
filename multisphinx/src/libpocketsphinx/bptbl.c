/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008-2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file bptbl.c Forward search lattice for N-Gram search.
 */

/* System headers. */
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "ngram_search.h"

bptbl_t *
bptbl_init(dict_t *dict, int n_alloc, int n_frame_alloc)
{
    bptbl_t *bptbl = ckd_calloc(1, sizeof(*bptbl));

    bptbl->dict = dict_retain(dict);
    bptbl->n_alloc = n_alloc;
    bptbl->n_frame_alloc = n_frame_alloc;

    bptbl->ent = ckd_calloc(bptbl->n_alloc, sizeof(*bptbl->ent));
    bptbl->word_idx = ckd_calloc(dict_size(dict),
                                 sizeof(*bptbl->word_idx));
    bptbl->bscore_stack_size = bptbl->n_alloc * 20;
    bptbl->bscore_stack = ckd_calloc(bptbl->bscore_stack_size,
                                     sizeof(*bptbl->bscore_stack));
    bptbl->ef_idx = ckd_calloc(bptbl->n_frame_alloc + 1,
                               sizeof(*bptbl->ef_idx));
    ++bptbl->ef_idx; /* Make bptableidx[-1] valid */
    bptbl->frm_wordlist = ckd_calloc(bptbl->n_frame_alloc,
                                     sizeof(*bptbl->frm_wordlist));

    return bptbl;
}

void
bptbl_free(bptbl_t *bptbl)
{
    if (bptbl == NULL)
        return;
    dict_free(bptbl->dict);
    ckd_free(bptbl->word_idx);
    ckd_free(bptbl->ent);
    ckd_free(bptbl->bscore_stack);
    if (bptbl->ef_idx != NULL)
        ckd_free(bptbl->ef_idx - 1);
    ckd_free(bptbl->frm_wordlist);
    ckd_free(bptbl);
}

void
dump_bptable(bptbl_t *bptbl)
{
    int i;
    E_INFO("Backpointer table (%d entries):\n", bptbl->n_ent);
    for (i = 0; i < bptbl->n_ent; ++i) {
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                    i, dict_wordstr(bptbl->dict, bptbl->ent[i].wid),
                    bptbl->ent[i].bp == -1 ? 0 : 
                    bptbl->ent[bptbl->ent[i].bp].frame + 1,
                    bptbl->ent[i].frame,
                    bptbl->ent[i].score,
                    bptbl->ent[i].bp);
    }
}

int
bptbl_push_frame(bptbl_t *bptbl, int oldest_bp, int frame_idx)
{
    E_INFO("pushing frame %d, oldest bp %d in frame %d\n",
           frame_idx, oldest_bp, oldest_bp == NO_BP ? -1 : bptbl->ent[oldest_bp].frame);
    if (oldest_bp >= 0)
        bptbl->window_sf = bptbl->ent[oldest_bp].frame + 1;
    if (frame_idx >= bptbl->n_frame_alloc) {
        bptbl->n_frame_alloc *= 2;
        bptbl->ef_idx = ckd_realloc(bptbl->ef_idx - 1,
                                    (bptbl->n_frame_alloc + 1)
                                    * sizeof(*bptbl->ef_idx));
        bptbl->frm_wordlist = ckd_realloc(bptbl->frm_wordlist,
                                          bptbl->n_frame_alloc
                                          * sizeof(*bptbl->frm_wordlist));
        ++bptbl->ef_idx; /* Make bptableidx[-1] valid */
    }
    bptbl->ef_idx[frame_idx] = bptbl->n_ent;
    return bptbl->n_ent;
}

void
bptable_gc(bptbl_t *bpt, int oldest_bp, int frame_idx)
{
    return;
#if 0
    E_INFO("Before sorting (%d : %d)\n", bpt->first_unsorted, oldest_bp);
    for (i = bpt->first_unsorted; i < oldest_bp; ++i) {
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                    i, dict_wordstr(bpt->dict, bpt->ent[i].wid),
                    bp_sf(bpt, i),
                    bpt->ent[i].frame,
                    bpt->ent[i].score,
                    bpt->ent[i].bp);
    }
    /* Insertion sort everything from first_unsorted to oldest_bp.
     * Although this is O(N**2), the backpointers are mostly sorted,
     * and we expect N to be fairly small */
    while (bpt->first_unsorted < (oldest_bp - 1)) {
        int i = bpt->first_unsorted, j;
        int sf;
        bp_t tmp;

        /* Short-cut already-sorted case. */
        sf = bp_sf(bpt, i);
        if (sf <= bp_sf(bpt, i + 1)) {
            ++bpt->first_unsorted;
            continue;
        }
        /* Copy this backpointer and it into the correct position. */
        memcpy(&tmp, bpt->ent + i, sizeof(tmp));
        while (i < (oldest_bp - 1) && sf > bp_sf(bpt, i + 1)) {
            memcpy(bpt->ent + i + 1, bpt->ent + i, sizeof(*bpt->ent));
            ++i;
        }
        memcpy(bpt->ent + i, &tmp, sizeof(tmp));
        /* Update all backpointers. */
        for (j = bpt->first_unsorted; j < oldest_bp; ++j) {
            if (bpt->ent[j].bp == bpt->first_unsorted)
                bpt->ent[j].bp = i;
            else if (bpt->ent[j].bp <= i)
                --bpt->ent[j].bp;
        }
        /* Increment first_unsorted. */
        ++bpt->first_unsorted;
    }
    bpt->first_unsorted = oldest_bp;
    E_INFO("After sorting (%d : %d)\n", old_first_unsorted, oldest_bp);
    for (i = old_first_unsorted; i < oldest_bp; ++i) {
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                    i, dict_wordstr(bpt->dict, bpt->ent[i].wid),
                    bp_sf(bpt, i),
                    bpt->ent[i].frame,
                    bpt->ent[i].score,
                    bpt->ent[i].bp);
    }
#endif
}
