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
bptbl_init(dict2pid_t *d2p, int n_alloc, int n_frame_alloc)
{
    bptbl_t *bptbl = ckd_calloc(1, sizeof(*bptbl));

    bptbl->d2p = dict2pid_retain(d2p);
    bptbl->n_alloc = n_alloc;
    bptbl->n_frame_alloc = n_frame_alloc;

    bptbl->ent = ckd_calloc(bptbl->n_alloc, sizeof(*bptbl->ent));
    bptbl->word_idx = ckd_calloc(dict_size(d2p->dict),
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
    dict2pid_free(bptbl->d2p);
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
                    i, dict_wordstr(bptbl->d2p->dict, bptbl->ent[i].wid),
                    bp_sf(bptbl, i),
                    bptbl->ent[i].frame,
                    bptbl->ent[i].score,
                    bptbl->ent[i].bp);
    }
}


#if 0
void
bptable_gc(bptbl_t *bpt, int oldest_bp, int frame_idx)
{
    E_INFO("Before sorting (%d : %d)\n", bpt->first_unsorted, oldest_bp);
    for (i = bpt->first_unsorted; i < oldest_bp; ++i) {
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                    i, dict_wordstr(bpt->d2p->dict, bpt->ent[i].wid),
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
                    i, dict_wordstr(bpt->d2p->dict, bpt->ent[i].wid),
                    bp_sf(bpt, i),
                    bpt->ent[i].frame,
                    bpt->ent[i].score,
                    bpt->ent[i].bp);
    }
}

static void
bptbl_gc(bptbl_t *bptbl, int oldest_bp)
{
    int prev_window_sf, window_sf;
    int32 *agenda;
    int i, j, n_bp, n_bp_active;

    /* window_sf is the first frame which is still active in search
     * (i.e. for which outgoing word arcs can still be generated).
     * Therefore, any future backpointer table entries will not point
     * backwards to any backpointers before (window_sf - 1), and thus
     * any backpointers which are not reachable from those exiting in
     * (window_sf - 1) will never be reachable. */
    prev_window_sf = bptbl->window_sf;
    window_sf = bptbl->ent[oldest_bp].frame + 1;
    E_INFO("window_sf %d prev_window_sf %d\n", window_sf, prev_window_sf);
    assert(window_sf >= prev_window_sf);
    if (window_sf == prev_window_sf)
        return;
    /* We can't garbage collect anything in window_sf - 1 since there
     * can still be backpointers pointing at it.  */
    if (window_sf - 1 == prev_window_sf)
        return;
    E_INFO("Garbage collecting from %d to %d\n",
           prev_window_sf, window_sf - 1);
    /* Invalidate all backpointer entries up to window_sf - 1. */
    for (i = 0; i < bptbl->ef_idx[window_sf - 1]; ++i)
        bptbl->ent[i].valid = FALSE;
    /* Now re-activate all reachable ones. */
    n_bp = bptbl->ef_idx[window_sf] - bptbl->ef_idx[window_sf - 1];
    E_INFO("Finding accessible from %d backpointers:\n", n_bp);
    /* The bptbl is a tree so there will never be more active nodes
     * than this in the agenda. */
    agenda = ckd_calloc(n_bp, sizeof(*agenda));
    for (i = 0; i < n_bp; ++i) {
        int32 bp = bptbl->ef_idx[bptbl->window_sf - 1] + i;
        agenda[i] = bp;
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
               bp, dict_wordstr(bptbl->d2p->dict,
                                bptbl->ent[bp].wid),
               bp_sf(bptbl, bp),
               bptbl->ent[bp].frame,
               bptbl->ent[bp].score,
               bptbl->ent[bp].bp);
    }
    n_bp_active = n_bp;
    while (n_bp_active) {
        E_INFO_NOFN("");
        for (j = 0; j < n_bp; ++j) {
            E_INFOCONT(" %-5d", agenda[j]);
            if (agenda[j] == NO_BP)
                continue;
            agenda[j] = bptbl->ent[agenda[j]].bp;
            if (agenda[j] == NO_BP)
                --n_bp_active;
            else {
                bptbl->ent[agenda[j]].valid = TRUE;
            }
        }
        E_INFOCONT("\n");
    }
    ckd_free(agenda);
    j = 0;
    for (i = 0; i < bptbl->ef_idx[window_sf - 1]; ++i)
        if (bptbl->ent[i].valid == FALSE)
            ++j;
    E_INFO("Invalidated %d entries\n", j);
    bptbl->window_sf = window_sf;
}
#endif

int
bptbl_push_frame(bptbl_t *bptbl, int oldest_bp, int frame_idx)
{
    E_INFO("pushing frame %d, oldest bp %d in frame %d\n",
           frame_idx, oldest_bp, oldest_bp == NO_BP ? -1 : bptbl->ent[oldest_bp].frame);
    bptbl->ef_idx[frame_idx] = bptbl->n_ent;
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
    return bptbl->n_ent;
}

bp_t *
bptbl_enter(bptbl_t *bptbl, int32 w, int frame_idx, int32 path,
            int32 score, int rc)
{
    int32 i, rcsize, *bss;
    bp_t *be;

    /* This might happen if recognition fails. */
    if (bptbl->n_ent == NO_BP) {
        E_ERROR("No entries in backpointer table!");
        return NULL;
    }

    /* Expand the backpointer tables if necessary. */
    if (bptbl->n_ent >= bptbl->n_alloc) {
        bptbl->n_alloc *= 2;
        bptbl->ent = ckd_realloc(bptbl->ent,
                                 bptbl->n_alloc
                                 * sizeof(*bptbl->ent));
        E_INFO("Resized backpointer table to %d entries\n", bptbl->n_alloc);
    }
    if (bptbl->bss_head >= bptbl->bscore_stack_size
        - bin_mdef_n_ciphone(bptbl->d2p->mdef)) {
        bptbl->bscore_stack_size *= 2;
        bptbl->bscore_stack = ckd_realloc(bptbl->bscore_stack,
                                          bptbl->bscore_stack_size
                                          * sizeof(*bptbl->bscore_stack));
        E_INFO("Resized score stack to %d entries\n", bptbl->bscore_stack_size);
    }

    bptbl->word_idx[w] = bptbl->n_ent;
    be = &(bptbl->ent[bptbl->n_ent]);
    be->wid = w;
    be->frame = frame_idx;
    be->bp = path;
    be->score = score;
    be->s_idx = bptbl->bss_head;
    be->valid = TRUE;
    bptbl_fake_lmstate(bptbl, bptbl->n_ent);

    /* DICT2PID */
    /* Get diphone ID for final phone and number of ssids corresponding to it. */
    be->last_phone = dict_last_phone(bptbl->d2p->dict,w);
    if (dict_is_single_phone(bptbl->d2p->dict, w)) {
        be->last2_phone = -1;
        rcsize = 1;
    }
    else {
        be->last2_phone = dict_second_last_phone(bptbl->d2p->dict, w);
        rcsize = dict2pid_rssid(bptbl->d2p, be->last_phone, be->last2_phone)->n_ssid;
    }
    /* Allocate some space on the bptbl->bscore_stack for all of these triphones. */
    for (i = rcsize, bss = bptbl->bscore_stack + bptbl->bss_head; i > 0; --i, bss++)
        *bss = WORST_SCORE;
    bptbl->bscore_stack[bptbl->bss_head + rc] = score;
    E_DEBUG(2,("Entered bp %d sf %d ef %d window_sf %d\n", bptbl->n_ent,
               bp_sf(bptbl, bptbl->n_ent), frame_idx, bptbl->window_sf));
    assert(bp_sf(bptbl, bptbl->n_ent) >= bptbl->window_sf);

    bptbl->n_ent++;
    bptbl->bss_head += rcsize;

    return be;
}

void
bptbl_fake_lmstate(bptbl_t *bptbl, int32 bp)
{
    int32 w, prev_bp;
    bp_t *be;

    assert(bp != NO_BP);

    be = &(bptbl->ent[bp]);
    prev_bp = bp;
    w = be->wid;

    while (dict_filler_word(bptbl->d2p->dict, w)) {
        prev_bp = bptbl->ent[prev_bp].bp;
        if (prev_bp == NO_BP)
            return;
        w = bptbl->ent[prev_bp].wid;
    }

    be->real_wid = dict_basewid(bptbl->d2p->dict, w);

    prev_bp = bptbl->ent[prev_bp].bp;
    be->prev_real_wid =
        (prev_bp != NO_BP) ? bptbl->ent[prev_bp].real_wid : -1;
}
