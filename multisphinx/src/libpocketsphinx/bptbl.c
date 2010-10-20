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
    bptbl->n_active_alloc = n_frame_alloc;

    bptbl->ent = ckd_calloc(bptbl->n_alloc, sizeof(*bptbl->ent));
    bptbl->permute = ckd_calloc(bptbl->n_alloc, sizeof(*bptbl->permute));
    bptbl->orig_sf = ckd_calloc(bptbl->n_alloc, sizeof(*bptbl->orig_sf));
    bptbl->word_idx = ckd_calloc(dict_size(d2p->dict),
                                 sizeof(*bptbl->word_idx));
    bptbl->bscore_stack_size = bptbl->n_alloc * 20;
    bptbl->bscore_stack = ckd_calloc(bptbl->bscore_stack_size,
                                     sizeof(*bptbl->bscore_stack));
    bptbl->ef_idx = ckd_calloc(bptbl->n_active_alloc,
                               sizeof(*bptbl->ef_idx));
    bptbl->frm_wordlist = ckd_calloc(bptbl->n_active_alloc,
                                     sizeof(*bptbl->frm_wordlist));
    bptbl->valid_fr = bitvec_alloc(bptbl->n_active_alloc);

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
    ckd_free(bptbl->permute);
    ckd_free(bptbl->orig_sf);
    ckd_free(bptbl->bscore_stack);
    ckd_free(bptbl->ef_idx);
    ckd_free(bptbl->frm_wordlist);
    bitvec_free(bptbl->valid_fr);
    ckd_free(bptbl);
}

void
bptbl_reset(bptbl_t *bptbl)
{
    int i;

    for (i = 0; i < bptbl->n_active_alloc; ++i) {
        bptbl->ef_idx[i] = -1;
    }
    bitvec_clear_all(bptbl->valid_fr, bptbl->n_active_alloc);
    bptbl->first_invert_bp = 0;
    bptbl->n_frame = 0;
    bptbl->n_ent = 0;
    bptbl->bss_head = 0;
    bptbl->active_fr = 0;
}

void
dump_bptable(bptbl_t *bptbl, int start, int end)
{
    int i, j;

    if (end == -1)
        end = bptbl->n_ent;
    E_INFO("Backpointer table (%d entries from %d to %d):\n",
           bptbl->n_ent, start, end);
    for (j = i = start; i < end; ++i) {
        if (bptbl->ent[i].valid == FALSE)
            continue;
        ++j;
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                    i, dict_wordstr(bptbl->d2p->dict, bptbl->ent[i].wid),
                    bp_sf(bptbl, i),
                    bptbl->ent[i].frame,
                    bptbl->ent[i].score,
                    bptbl->ent[i].bp);
    }
    E_INFO("%d valid entries\n", j);
    if (bptbl_ef_idx(bptbl, start) < bptbl->active_fr)
        start = bptbl_ef_idx(bptbl, bptbl->active_fr);
    E_INFO("End frame index from %d to %d:\n",
           bptbl->ent[start].frame, bptbl->ent[end-1].frame);
    for (i = bptbl->ent[start].frame; i <= bptbl->ent[end-1].frame; ++i) {
        E_INFO_NOFN("%d: %d\n", i, bptbl_ef_idx(bptbl, i));
        assert(bptbl_ef_idx(bptbl, i) >= bptbl_ef_idx(bptbl, i-1));
    }
}

/**
 * Mark coaccessible entries in the backpointer table.
 *
 * @param bptbl Backpointer table
 * @param sf First frame (as per ef_idx) to mark entries in
 * @param ef Frame after last frame (as per ef_idx) to mark entries in
 * @param cf Current frame of search
 */
static void
bptbl_mark(bptbl_t *bptbl, int sf, int ef, int cf)
{
    int i, j, n_active_fr, last_gc_fr;

    /* Invalidate all backpointer entries between sf and ef. */
    E_DEBUG(2,("Garbage collecting from %d to %d (%d to %d):\n",
               bptbl_ef_idx(bptbl, sf), bptbl_ef_idx(bptbl, ef), sf, ef));
    for (i = bptbl_ef_idx(bptbl, sf);
         i < bptbl_ef_idx(bptbl, ef); ++i) {
        bptbl->ent[i].valid = FALSE;
    }

    /* Now re-activate all ones backwards reachable from the elastic
     * window. (make sure cf has been pushed!) */
    E_DEBUG(2, ("Finding coaccessible frames from backpointers from %d to %d\n",
                bptbl_ef_idx(bptbl, ef), bptbl_ef_idx(bptbl, cf)));
    /* Collect coaccessible frames from these backpointers */
    bitvec_clear_all(bptbl->valid_fr, cf - bptbl->active_fr);
    n_active_fr = 0;
    for (i = bptbl_ef_idx(bptbl, ef);
         i < bptbl_ef_idx(bptbl, cf); ++i) {
        int32 bp = bptbl->ent[i].bp;
        if (bp != NO_BP) {
            int frame = bptbl->ent[bp].frame;
            if (frame >= bptbl->active_fr
                && bitvec_is_clear(bptbl->valid_fr, frame - bptbl->active_fr)) {
                bitvec_set(bptbl->valid_fr, frame - bptbl->active_fr);
                ++n_active_fr;
            }
        }
    }
    /* Track the last frame with outgoing backpointers for gc */
    last_gc_fr = ef - 1;
    while (n_active_fr > 0) {
        int next_gc_fr = 0;
        n_active_fr = 0;
        for (i = sf; i <= last_gc_fr; ++i) {
            if (i >= bptbl->active_fr
                && bitvec_is_set(bptbl->valid_fr, i - bptbl->active_fr)) {
                bitvec_clear(bptbl->valid_fr, i - bptbl->active_fr);
                /* Add all backpointers in this frame (the bogus
                 * lattice generation algorithm) */
                for (j = bptbl_ef_idx(bptbl, i);
                     j < bptbl_ef_idx(bptbl, i + 1); ++j) {
                    int32 bp = bptbl->ent[j].bp;
                    bptbl->ent[j].valid = TRUE;
                    if (bp != NO_BP) {
                        int frame = bptbl->ent[bp].frame;
                        if (frame >= bptbl->active_fr
                            && bitvec_is_clear(bptbl->valid_fr,
                                               frame - bptbl->active_fr)) {
                            bitvec_set(bptbl->valid_fr,
                                       frame - bptbl->active_fr);
                            ++n_active_fr;
                        }
                        if (frame > next_gc_fr)
                            next_gc_fr = frame;
                    }
                }
            }
        }
        E_DEBUG(3,("last_gc_fr %d => %d\n", last_gc_fr, next_gc_fr));
        last_gc_fr = next_gc_fr;
    }
    E_DEBUG(2,("Removed"));
    for (i = bptbl_ef_idx(bptbl, sf);
         i < bptbl_ef_idx(bptbl, ef); ++i) {
        if (bptbl->ent[i].valid == FALSE)
            E_DEBUGCONT(2,(" %d", i));
    }
    E_DEBUGCONT(2,("\n"));
}

/**
 * Compact the backpointer table.
 *
 * @param bptbl Backpointer table.
 * @param eidx First index which cannot be compacted.
 * @return End of compacted backpointers.
 */
static int
bptbl_compact(bptbl_t *bptbl, int eidx)
{
    int src, dest, ef;
    int sidx = bptbl->first_invert_bp;

    ef = bptbl->ent[sidx].frame;
    E_INFO("compacting from %d to %d (%d to %d)\n",
           sidx, eidx, ef, bptbl->ent[eidx].frame);
    for (dest = src = sidx; src < eidx; ++src) {
        /* Update all ef_idx including missing frames (there are many) */
        while (ef >= bptbl->active_fr && ef <= bptbl->ent[src].frame) {
            assert(ef - bptbl->active_fr < bptbl->n_active_alloc);
            bptbl->ef_idx[ef - bptbl->active_fr] = dest;
            ++ef;
        }
        if (bptbl->ent[src].valid) {
            if (dest != src) {
                bptbl->ent[dest] = bptbl->ent[src];
                bptbl->ent[src].valid = FALSE;
            }
            E_DEBUG(4,("permute %d => %d\n", src, dest));
            bptbl->permute[src] = dest;
            ++dest;
        }
        else {
            E_DEBUG(4,("permute %d => -1\n", src));
            bptbl->permute[src] = -1;
        }
    }
    if (eidx == bptbl->n_ent)
        bptbl->n_ent = dest;
    return dest;
}

/**
 * Sort retired backpointer entries by start frame.
 *
 * @param bptbl Backpointer table.
 * @param sidx Index of first unordered entry.
 * @param eidx Index of last unordered entry plus one.
 */
static int
bptbl_forward_sort(bptbl_t *bptbl, int sidx, int eidx)
{
    int i;

    E_INFO("Sorting forward from %d to %d\n", sidx, eidx);
    /* Straightforward for now, we just insertion sort these dudes and
     * update the permutation table and first_invert_bp as necessary.
     * This could be done in conjunction with compaction but it
     * becomes very complicated and possibly slower. */

    /* Reset the permutation table (for now we have to do the
     * backpointer update twice, we'll fix it later by storing the
     * reverse permutation table and explicitly inverting it) */
    for (i = 0; i < bptbl->n_ent; ++i) {
        bptbl->permute[i] = i;
        bptbl->orig_sf[i] = bp_sf(bptbl, i);
    }

    /* Crawl from sidx to eidx inserting and updating the permutation
     * table as necessary. */
    for (i = sidx; i < eidx; ++i) {
        bp_t ent;
        int j, k, isf;
        isf = bptbl->orig_sf[i];
        for (j = i - 1; j >= 0; --j)
            if (bptbl->orig_sf[j] <= isf)
                break;
        ++j;
        if (j == i)
            continue;
        E_INFO("Inserting %d (sf %d) to %d (sf %d)\n",
               i, isf, j, bptbl->orig_sf[j]);
        ent = bptbl->ent[i];
        memmove(bptbl->ent + j + 1, bptbl->ent + j,
                (i - j) * sizeof(*bptbl->ent));
        memmove(bptbl->orig_sf + j + 1, bptbl->orig_sf + j,
                (i - j) * sizeof(*bptbl->orig_sf));
        bptbl->ent[j] = ent;
        bptbl->orig_sf[j] = isf;
        bptbl->permute[i] = j;
        E_INFO("permute %d => %d\n", i, j);
        for (k = 0; k < i; ++k) {
            if (bptbl->permute[k] >= j && bptbl->permute[k] < i) {
                E_INFO("permute %d => %d\n", k, bptbl->permute[k]+1);
                ++bptbl->permute[k];
            }
        }
    }

    return 0;
}

/**
 * Remap backpointers in backpointer table.
 *
 * @param bptbl Backpointer table.
 * @param sidx Index of first backpointer to be updated
 * @param eidx Index of last backpointer to be updated plus one.
 * @param last_invert_bp Index of first backpointer which can be remapped plus one.
 * @param last_invert_bp Index of last backpointer which can be remapped plus one.
 */
static void
bptbl_invert(bptbl_t *bptbl, int sidx, int eidx,
             int first_invert_bp, int last_invert_bp)
{
    int i;

    E_INFO("inverting %d:%d from %d to %d\n", first_invert_bp, last_invert_bp, sidx, eidx);
    for (i = sidx; i < eidx; ++i) {
        if (bptbl->ent[i].bp >= first_invert_bp && bptbl->ent[i].bp < last_invert_bp) {
            if (bptbl->ent[i].bp != bptbl->permute[bptbl->ent[i].bp])
                E_INFO("invert %d => %d in %d\n",
                       bptbl->ent[i].bp, bptbl->permute[bptbl->ent[i].bp], i);
            bptbl->ent[i].bp = bptbl->permute[bptbl->ent[i].bp];
            assert(bp_sf(bptbl, i) < bptbl->ent[i].frame);
        }
    }
}

/**
 * Update the active frame pointer.
 */
static void
bptbl_update_active_fr(bptbl_t *bptbl, int active_fr)
{
    int delta = active_fr - bptbl->active_fr;
    if (delta == 0)
        return;
    assert(delta > 0);
    E_DEBUG(3,("moving %d ef_idx from %d (%d - %d)\n",
               bptbl->n_frame - active_fr,
               delta, active_fr, bptbl->active_fr));
    assert(delta + bptbl->n_frame - active_fr < bptbl->n_active_alloc);
    memmove(bptbl->ef_idx, bptbl->ef_idx + delta,
            (bptbl->n_frame - active_fr) * sizeof(*bptbl->ef_idx));
    bptbl->active_fr = active_fr;
}

static void
bptbl_gc(bptbl_t *bptbl, int oldest_bp, int frame_idx)
{
    int prev_active_fr, active_fr;
    int last_compacted_bp;

    /* active_fr is the first frame which is still active in search
     * (i.e. for which outgoing word arcs can still be generated).
     * Therefore, any future backpointer table entries will not point
     * backwards to any backpointers before active_fr, and thus any
     * backpointers which are not reachable from those exiting in
     * active_fr will never be reachable. */
    prev_active_fr = bptbl->active_fr;
    active_fr = (oldest_bp == NO_BP) ? 0 : bptbl->ent[oldest_bp].frame;
    assert(active_fr >= prev_active_fr);
    if (active_fr <= prev_active_fr + 1)
        return;
    /* If there is nothing to GC then finish up */
    if (bptbl_ef_idx(bptbl, prev_active_fr) == bptbl_ef_idx(bptbl, active_fr)) {
        bptbl_update_active_fr(bptbl, active_fr);
        return;
    }

    bptbl_mark(bptbl, prev_active_fr, active_fr, frame_idx);
#if 0
    E_INFO("before compaction\n");
    dump_bptable(bptbl, 0, -1);
#endif
    last_compacted_bp = bptbl_compact(bptbl, bptbl_ef_idx(bptbl, active_fr));

    E_INFO("after compaction\n");
    dump_bptable(bptbl, 0, -1);
    bptbl_invert(bptbl, bptbl->first_invert_bp, last_compacted_bp,
                 bptbl->first_invert_bp, bptbl_ef_idx(bptbl, active_fr));
    bptbl_invert(bptbl, bptbl_ef_idx(bptbl, active_fr), bptbl->n_ent,
                 bptbl->first_invert_bp, bptbl_ef_idx(bptbl, active_fr));
    E_INFO("after inversion\n");
    dump_bptable(bptbl, 0, -1);
    bptbl_forward_sort(bptbl, bptbl->first_invert_bp, last_compacted_bp);
    E_INFO("after sort\n");
    dump_bptable(bptbl, 0, -1);
    bptbl_invert(bptbl, bptbl->first_invert_bp, last_compacted_bp,
                 /* FIXME: don't actually have to start at 0. */
                 0, bptbl_ef_idx(bptbl, active_fr));
    bptbl_invert(bptbl, bptbl_ef_idx(bptbl, active_fr), bptbl->n_ent,
                 /* FIXME: don't actually have to start at 0. */
                 0, bptbl_ef_idx(bptbl, active_fr));
    E_INFO("after inversion 2\n");
    dump_bptable(bptbl, 0, -1);
    int i;
    for (i = bptbl->first_invert_bp; i < last_compacted_bp; ++i) {
        E_INFO_NOFN("%d: orig_sf %d bp_sf %d\n",
                    i, bptbl->orig_sf[i], bp_sf(bptbl, i));
        assert(bptbl->orig_sf[i] == bp_sf(bptbl, i));
    }
    for (i = bptbl_ef_idx(bptbl, active_fr); i < bptbl->n_ent; ++i) {
        E_INFO_NOFN("%d: orig_sf %d bp_sf %d\n",
                    i, bptbl->orig_sf[i], bp_sf(bptbl, i));
        assert(bptbl->orig_sf[i] == bp_sf(bptbl, i));
    }
    bptbl->first_invert_bp = last_compacted_bp;
    bptbl_update_active_fr(bptbl, active_fr);
}

int
bptbl_push_frame(bptbl_t *bptbl, int oldest_bp, int frame_idx)
{
    E_DEBUG(2,("pushing frame %d, oldest bp %d in frame %d\n",
               frame_idx, oldest_bp, oldest_bp == NO_BP ? -1 : bptbl->ent[oldest_bp].frame));
    if (frame_idx - bptbl->active_fr >= bptbl->n_active_alloc) {
        bptbl->n_active_alloc *= 2;
        E_INFO("Reallocating frame-based bptr arrays to %d\n", bptbl->n_active_alloc);
        bptbl->ef_idx = ckd_realloc(bptbl->ef_idx,
                                    bptbl->n_active_alloc * sizeof(*bptbl->ef_idx));
        bptbl->frm_wordlist = ckd_realloc(bptbl->frm_wordlist,
                                          bptbl->n_active_alloc
                                          * sizeof(*bptbl->frm_wordlist));
        bptbl->valid_fr = bitvec_realloc(bptbl->valid_fr, bptbl->n_active_alloc);
    }
    bptbl->ef_idx[frame_idx - bptbl->active_fr] = bptbl->n_ent;
    bptbl->n_frame = frame_idx + 1;
    bptbl_gc(bptbl, oldest_bp, frame_idx);
    return bptbl->n_ent;
}

int32
bptbl_ef_idx(bptbl_t *bptbl, int frame_idx)
{
    if (frame_idx < bptbl->active_fr)
        return 0;
    else if (frame_idx >= bptbl->n_frame)
        return bptbl->n_ent;
    else {
        return bptbl->ef_idx[frame_idx - bptbl->active_fr];
    }
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
        bptbl->permute = ckd_realloc(bptbl->permute,
                                     bptbl->n_alloc
                                     * sizeof(*bptbl->permute));
        bptbl->orig_sf = ckd_realloc(bptbl->orig_sf,
                                     bptbl->n_alloc
                                     * sizeof(*bptbl->orig_sf));
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
    E_DEBUG(3,("Entered bp %d sf %d ef %d active_fr %d\n", bptbl->n_ent,
               bp_sf(bptbl, bptbl->n_ent), frame_idx, bptbl->active_fr));
    assert(bp_sf(bptbl, bptbl->n_ent) >= bptbl->active_fr);

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
