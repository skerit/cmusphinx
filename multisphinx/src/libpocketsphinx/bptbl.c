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
#include "bptbl.h"

/* Special unlocked versions of some utility functions */
static int bptbl_rcsize(bptbl_t *bptbl, bp_t *be);
static bpidx_t bptbl_find_exit(bptbl_t *bptbl, int32 wid);
static int32 bptbl_ef_idx_internal(bptbl_t *bptbl, int frame_idx);
static bp_t *bptbl_ent_internal(bptbl_t *bptbl, bpidx_t bpidx);
static int bptbl_sf_internal(bptbl_t *bptbl, bpidx_t bpidx);
static int bptbl_ef_count_internal(bptbl_t *bptbl, int frame_idx);

#if 0
#undef E_DEBUG
#define E_DEBUG(level,x) E_INFO x
#undef E_DEBUGCONT
#define E_DEBUGCONT(level,x) E_INFOCONT x
#endif

bptbl_t *
bptbl_init(dict2pid_t *d2p, int n_alloc, int n_frame_alloc)
{
    bptbl_t *bptbl = ckd_calloc(1, sizeof(*bptbl));

    bptbl->d2p = dict2pid_retain(d2p);

    bptbl->ent = garray_init(0, sizeof(bp_t));
    garray_reserve(bptbl->ent, n_alloc / 2);
    bptbl->retired = garray_init(0, sizeof(bp_t));
    garray_reserve(bptbl->retired, n_alloc / 2);
    bptbl->permute = garray_init(0, sizeof(bpidx_t));
    garray_reserve(bptbl->permute, n_frame_alloc);
    bptbl->ef_idx = garray_init(0, sizeof(bpidx_t));
    garray_reserve(bptbl->ef_idx, n_frame_alloc);
    bptbl->n_frame_alloc = n_frame_alloc;
    bptbl->valid_fr = bitvec_alloc(bptbl->n_frame_alloc);
    bptbl->rc = garray_init(0, sizeof(int32));
    garray_reserve(bptbl->rc, n_alloc * 20); /* 20 = guess at average number of rcs/word */

    bptbl->evt = sbevent_init(FALSE);
    bptbl->mtx = sbmtx_init();
    return bptbl;
}

bptbl_t *
bptbl_retain(bptbl_t *bpt)
{
    ++bpt->refcount;
    return bpt;
}

int
bptbl_free(bptbl_t *bptbl)
{
    if (bptbl == NULL)
        return 0;
    if (--bptbl->refcount > 0)
        return bptbl->refcount;

    dict2pid_free(bptbl->d2p);
    garray_free(bptbl->ent);
    garray_free(bptbl->retired);
    garray_free(bptbl->permute);
    garray_free(bptbl->ef_idx);
    garray_free(bptbl->rc);
    bitvec_free(bptbl->valid_fr);
    sbevent_free(bptbl->evt);
    sbmtx_free(bptbl->mtx);
    ckd_free(bptbl);
    return 0;
}

void
bptbl_reset(bptbl_t *bptbl)
{
    sbmtx_lock(bptbl->mtx);
    bitvec_clear_all(bptbl->valid_fr, bptbl->n_frame_alloc);
    garray_reset(bptbl->ent);
    garray_reset(bptbl->permute);
    garray_reset(bptbl->ef_idx);
    garray_reset(bptbl->retired);
    garray_reset(bptbl->rc);
    bptbl->dest_s_idx = 0;
    bptbl->n_frame = 0;
    bptbl->oldest_bp = NO_BP;
    sbmtx_unlock(bptbl->mtx);
}

void
bptbl_dump(bptbl_t *bptbl)
{
    int i;

    sbmtx_lock(bptbl->mtx);
    E_INFO("Retired backpointers (%d entries, oldest active %d):\n",
           bptbl_retired_idx(bptbl), bptbl->oldest_bp);
    for (i = garray_base(bptbl->retired);
         i < bptbl_retired_idx(bptbl); ++i) {
        bp_t *ent = garray_ptr(bptbl->retired, bp_t, i);
        assert(ent->valid);
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                    i, dict_wordstr(bptbl->d2p->dict, ent->wid),
                    bptbl_sf_internal(bptbl, i),
                    ent->frame,
                    ent->score,
                    ent->bp);
    }
    E_INFO("Active backpointers (%d entries starting at %d):\n",
           bptbl_end_idx(bptbl) - bptbl_active_idx(bptbl),
           bptbl_active_idx(bptbl));
    for (i = bptbl_active_idx(bptbl); i < bptbl_end_idx(bptbl); ++i) {
        bp_t *ent = garray_ptr(bptbl->ent, bp_t, i);
        if (!ent->valid)
            E_INFO_NOFN("%-5d INVALID\n", i);
        else
            E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                        i, dict_wordstr(bptbl->d2p->dict, ent->wid),
                        bptbl_sf_internal(bptbl, i),
                        ent->frame,
                        ent->score,
                        ent->bp);
    }
    sbmtx_unlock(bptbl->mtx);
}

/**
 * Mark coaccessible active entries in the backpointer table.
 *
 * @param bptbl Backpointer table
 * @param ef Frame after last frame (as per ef_idx) to mark entries in
 * @param cf Current frame of search
 * @return Number of remaining active entries in marked region
 */
static int
bptbl_mark(bptbl_t *bptbl, int ef, int cf)
{
    int i, j, n_active_fr, last_gc_fr;

    assert(cf >= ef);
    assert(ef > bptbl_active_frame(bptbl));

    /* Invalidate all active backpointer entries up to ef. */
    E_DEBUG(2,("Invalidating backpointers from %d to %d (%d to %d)\n",
               bptbl_active_frame(bptbl), ef,
               bptbl_ef_idx_internal(bptbl, bptbl_active_frame(bptbl)),
               bptbl_ef_idx_internal(bptbl, ef)));
    for (i = bptbl_ef_idx_internal(bptbl, bptbl_active_frame(bptbl));
         i < bptbl_ef_idx_internal(bptbl, ef); ++i) {
        E_DEBUG(5,("Invalidate bp %d\n", i));
        garray_ent(bptbl->ent, bp_t, i).valid = FALSE;
    }

    /* Now re-activate all ones backwards reachable from the search graph. */
    E_DEBUG(2,("Finding coaccessible frames from backpointers from %d to %d (%d to %d)\n",
               ef, cf,
               bptbl_ef_idx_internal(bptbl, ef), bptbl_ef_idx_internal(bptbl, cf)));
    /* Mark everything immediately reachable from (ef..cf) */
    bitvec_clear_all(bptbl->valid_fr, cf - bptbl_active_frame(bptbl));
    n_active_fr = 0;
    /* NOTE: This for statement can be sped up at the cost of being
     * less obvious. */
    for (i = bptbl_ef_idx_internal(bptbl, ef);
         i < bptbl_ef_idx_internal(bptbl, cf); ++i) {
        bp_t *ent, *prev;
        ent = garray_ptr(bptbl->ent, bp_t, i);
        prev = garray_ptr(bptbl->ent, bp_t, ent->bp);
        if (!ent->valid) /* May be invalidated by maxwpf */
            continue;
        if (prev != NULL) {
            int frame = prev->frame;
            if (frame >= bptbl_active_frame(bptbl)
                && bitvec_is_clear(bptbl->valid_fr,
                                   frame - bptbl_active_frame(bptbl))) {
                E_DEBUG(5,("Validate frame %d\n", frame - bptbl_active_frame(bptbl)));
                bitvec_set(bptbl->valid_fr, frame - bptbl_active_frame(bptbl));
                ++n_active_fr;
            }
        }
    }
    
    /* Track the last frame with outgoing backpointers for gc */
    last_gc_fr = ef - 1;
    /* Walk back from every frame marked active up to the last one
     * marked active in the previous round (this means we scan some of
     * them multiple times, if this gets slow this can be
     * optimized). */
    while (n_active_fr > 0) {
        int next_gc_fr = 0;
        n_active_fr = 0;
        for (i = bptbl_active_frame(bptbl); i <= last_gc_fr; ++i) {
            /* NOTE: We asserted i >= bptbl_active_frame(bptbl) */
            if (bitvec_is_set(bptbl->valid_fr, i - bptbl_active_frame(bptbl))) {
                bitvec_clear(bptbl->valid_fr, i - bptbl_active_frame(bptbl));
                /* Add all backpointers in this frame (the bogus
                 * lattice generation algorithm) */
                for (j = bptbl_ef_idx_internal(bptbl, i);
                     j < bptbl_ef_idx_internal(bptbl, i + 1); ++j) {
                    bp_t *ent = garray_ptr(bptbl->ent, bp_t, j);
                    bp_t *prev = garray_ptr(bptbl->ent, bp_t, ent->bp);
                    ent->valid = TRUE;
                    if (prev != NULL) {
                        int frame = prev->frame;
                        if (frame >= bptbl_active_frame(bptbl)
                            && bitvec_is_clear(bptbl->valid_fr,
                                               frame - bptbl_active_frame(bptbl))) {
                            bitvec_set(bptbl->valid_fr,
                                       frame - bptbl_active_frame(bptbl));
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
    for (j = 0, i = bptbl_ef_idx_internal(bptbl, bptbl_active_frame(bptbl));
         i < bptbl_ef_idx_internal(bptbl, ef); ++i) {
        bp_t *ent = garray_ptr(bptbl->ent, bp_t, i);
        assert(ent != NULL);
        if (ent->valid)
            ++j;
        else
            E_DEBUGCONT(2,(" %d", i));
    }
    E_DEBUGCONT(2,("\n"));
    return j;
}

/**
 * Retire accessible backpointers
 *
 * @param bptbl Backpointer table.
 * @param n_bp Number of backpointers to be retired.
 * @param eidx Absolute index of first backpointer which cannot be retired.
 * @return Index of last retired backpointer (relative to bptbl->retired) plus one.
 */
static int
bptbl_retire(bptbl_t *bptbl, int n_retired, int eidx)
{
    int src, dest;
    /* bss index for active frames (we don't want to track this in bptbl). */
    int active_dest_s_idx;

    /* First available backpointer index in retired. */
    dest = bptbl_retired_idx(bptbl);
    /* Expand retired if necessary. */
    garray_expand_to(bptbl->retired, dest + n_retired);
    /* Note we use the "cooked" backpointer indices here. */
    for (src = bptbl_active_idx(bptbl); src < eidx; ++src) {
        bp_t *src_ent = garray_ptr(bptbl->ent, bp_t, src);
        bp_t *dest_ent = garray_ptr(bptbl->retired, bp_t, dest);
        if (src_ent->valid) {
            int rcsize = bptbl_rcsize(bptbl, src_ent);
            E_DEBUG(4,("permute %d => %d\n", src, dest));
            if (src_ent->s_idx != bptbl->dest_s_idx) {
                E_DEBUG(4,("Moving %d rc scores from %d to %d for bptr %d\n",
                           rcsize, src_ent->s_idx, bptbl->dest_s_idx, dest));
                assert(src_ent->s_idx > bptbl->dest_s_idx);
                if (src < bptbl_end_idx(bptbl) - 1) {
                    bp_t *src1_ent = garray_ptr(bptbl->ent, bp_t, src + 1);
                    assert(bptbl->dest_s_idx + rcsize <= src1_ent->s_idx);
                }
                garray_move(bptbl->rc, bptbl->dest_s_idx, src_ent->s_idx, rcsize);
            }
            *dest_ent = *src_ent;
            dest_ent->s_idx = bptbl->dest_s_idx;

            assert(src < garray_next_idx(bptbl->permute));
            garray_ent(bptbl->permute, bpidx_t, src) = dest;
            bptbl->dest_s_idx += rcsize;
            ++dest;
        }
        else {
            E_DEBUG(4,("permute %d => -1 src %d ef_idx %d\n",
                       src, src - bptbl_active_idx(bptbl), bptbl_active_idx(bptbl)));
            assert(src < garray_next_idx(bptbl->permute));
            garray_ent(bptbl->permute, bpidx_t, src) = -1;
        }
    }
    /* We can keep compacting the bscore_stack since it is indirected. */
    if (src < bptbl_end_idx(bptbl)
        && garray_ent(bptbl->ent, bp_t, src).s_idx != bptbl->dest_s_idx) {
        /* Leave dest_s_idx where it is for future compaction. */
        active_dest_s_idx = bptbl->dest_s_idx;
        while (src < bptbl_end_idx(bptbl)) {
            bp_t *src_ent = garray_ptr(bptbl->ent, bp_t, src);
            int rcsize = bptbl_rcsize(bptbl, src_ent);
            E_DEBUG(4,("Moving %d rc scores from %d to %d for bptr %d\n",
                       rcsize, src_ent->s_idx, active_dest_s_idx, src));
            if (src < bptbl_end_idx(bptbl) - 1) {
                bp_t *src1_ent = garray_ptr(bptbl->ent, bp_t, src + 1);
                assert(bptbl->dest_s_idx + rcsize <= src1_ent->s_idx);
            }
            garray_move(bptbl->rc, active_dest_s_idx, src_ent->s_idx, rcsize);
            src_ent->s_idx = active_dest_s_idx;
            active_dest_s_idx += rcsize;
            ++src;
        }
        garray_pop_from(bptbl->rc, active_dest_s_idx);
    }
    sbmtx_unlock(bptbl->mtx);
    return dest;
}


/**
 * Remap backpointer indices in backpointer table.
 *
 * This needs to be done both in newly retired backpointers (which
 * start at first_retired_bp) as well as in still active backpointers,
 * which range from first_active_bp to bptbl->n_ent.
 *
 * This has the side effect of setting bptbl->oldest_bp
 *
 * @param bptbl Backpointer table.
 * @param last_retired_bp Index of last retired backpointer plus one
 * @param last_remapped_bp Last backpointer index to be remapped
 * @param first_active_bp Index of first remaining active backpointer
 */
static void
bptbl_remap(bptbl_t *bptbl, int first_retired_bp,
            int last_remapped_bp, int first_active_bp)
{
    int last_retired_bp;
    int i;

    /* Need to protect bptbl->retired and bptbl->rc */
    sbmtx_lock(bptbl->mtx);
    last_retired_bp = bptbl_retired_idx(bptbl);
    E_DEBUG(2,("remapping %d:%d from %d to %d and %d to %d\n",
               first_retired_bp, last_remapped_bp,
               first_retired_bp, last_retired_bp,
               first_active_bp, bptbl_end_idx(bptbl)));
    /* First remap backpointers in newly retired bps. */
    for (i = first_retired_bp; i < last_retired_bp; ++i) {
        bp_t *bpe = garray_ptr(bptbl->retired, bp_t, i);
        /* Remember, these are the *source* backpointer indices, so
         * they fall in the range between prev_active_fr (which is the
         * first index of ef_idx) and active_fr. */
        if (bpe->bp >= bptbl_active_idx(bptbl)
            && bpe->bp < last_remapped_bp) {
            assert(bpe->bp < garray_next_idx(bptbl->permute));
            if (bpe->bp != garray_ent(bptbl->permute, bpidx_t, bpe->bp))
                E_DEBUG(4,("remap retired %d => %d in %d\n",
                           bpe->bp, garray_ent(bptbl->permute, bpidx_t, bpe->bp), i));
            bpe->bp = garray_ent(bptbl->permute, bpidx_t, bpe->bp);
            assert(bptbl_sf_internal(bptbl, i) <= bpe->frame);
        }
    }

    /* Now remap backpointers in still-active bps (which point to the
     * newly retired ones) */
    bptbl->oldest_bp = last_retired_bp - 1;
    for (i = first_active_bp; i < bptbl_end_idx(bptbl); ++i) {
        bp_t *bpe = garray_ptr(bptbl->ent, bp_t, i);
        if (bpe->bp >= bptbl_active_idx(bptbl) && bpe->bp < last_remapped_bp) {
            assert(bpe->bp < garray_next_idx(bptbl->permute));
            if (bpe->bp != garray_ent(bptbl->permute, bpidx_t, bpe->bp))
                E_DEBUG(4,("remap active %d => %d in %d\n",
                           bpe->bp, garray_ent(bptbl->permute, bpidx_t, bpe->bp), i));
            bpe->bp = garray_ent(bptbl->permute, bpidx_t, bpe->bp);
            assert(bptbl_sf_internal(bptbl, i) <= bpe->frame);
        }
        if (bpe->bp < bptbl->oldest_bp)
            bptbl->oldest_bp = bpe->bp;
    }
    sbmtx_unlock(bptbl->mtx);
}

/**
 * Update the active frame pointer and backpointer array.
 */
static void
bptbl_update_active(bptbl_t *bptbl, int active_fr)
{
    /* This means nothing happened. */
    if (active_fr == bptbl_active_frame(bptbl))
        return;

    /* Shift back active backpointers. */
    E_DEBUG(3,("moving %d ent from %d (%d - %d)\n",
               bptbl_end_idx(bptbl) - garray_ent(bptbl->ef_idx, bpidx_t, active_fr),
               garray_ent(bptbl->ef_idx, bpidx_t, active_fr) - bptbl_active_idx(bptbl),
               garray_ent(bptbl->ef_idx, bpidx_t, active_fr),
               bptbl_active_idx(bptbl)));
    garray_shift_from(bptbl->ent, garray_ent(bptbl->ef_idx, bpidx_t, active_fr));
    garray_set_base(bptbl->ent, garray_ent(bptbl->ef_idx, bpidx_t, active_fr));

    /* Shift back end frame indices (implicitly updating output of bptbl_active_idx) */
    E_DEBUG(3,("moving %d ef_idx from %d (%d - %d)\n",
               bptbl->n_frame - active_fr,
               active_fr - bptbl_active_frame(bptbl),
               active_fr, bptbl_active_frame(bptbl)));
    garray_shift_from(bptbl->ef_idx, active_fr);
    garray_set_base(bptbl->ef_idx, active_fr);
}

static void
bptbl_gc(bptbl_t *bptbl, int oldest_bp, int frame_idx)
{
    int next_active_fr;
    int n_retired;
    int first_retired_bp;

    /* active_fr is the first frame which is still active in search
     * (i.e. for which outgoing word arcs can still be generated).
     * Therefore, any future backpointer table entries will not point
     * backwards to any backpointers before active_fr, and thus any
     * backpointers which are not reachable from those exiting in
     * active_fr will never be reachable. */
    if (oldest_bp == NO_BP)
        next_active_fr = 0;
    else
        next_active_fr = garray_ent(bptbl->ent, bp_t, oldest_bp).frame;
    assert(next_active_fr >= bptbl_active_frame(bptbl));
    /* Need at least 2 frames to GC. */
    if (next_active_fr <= bptbl_active_frame(bptbl) + 1)
        return;
    /* If there is nothing to GC then finish up. */
    if (bptbl_ef_idx_internal(bptbl, bptbl_active_frame(bptbl))
        == bptbl_ef_idx_internal(bptbl, next_active_fr)) {
        bptbl_update_active(bptbl, next_active_fr);
        return;
    }
    E_DEBUG(2,("GC from frame %d to %d\n", bptbl_active_frame(bptbl),
               next_active_fr));
    /* Expand the permutation table if necessary. */
    garray_expand_to(bptbl->permute, bptbl_ef_idx_internal(bptbl, next_active_fr));
    garray_set_base(bptbl->permute, bptbl_active_idx(bptbl));
    /* Mark, compact, snap pointers. */
    n_retired = bptbl_mark(bptbl, next_active_fr, frame_idx);
    E_DEBUG(2,("About to retire %d bps\n", n_retired));
    first_retired_bp = bptbl_retired_idx(bptbl);
    bptbl_retire(bptbl, n_retired,
                 bptbl_ef_idx_internal(bptbl, next_active_fr));
    bptbl_remap(bptbl, first_retired_bp,
                bptbl_ef_idx_internal(bptbl, next_active_fr),
                bptbl_ef_idx_internal(bptbl, next_active_fr));
    bptbl_update_active(bptbl, next_active_fr);
    E_DEBUG(2,("Retired %d bps: now %d retired, %d active\n", n_retired,
               bptbl_retired_idx(bptbl),
               bptbl_end_idx(bptbl) - bptbl_active_idx(bptbl)));
    E_DEBUG(2,("First active sf %d output frame %d window %d\n",
           bptbl->oldest_bp == NO_BP
           ? 0 : garray_ent(bptbl->retired, bp_t, bptbl->oldest_bp).frame + 1,
           frame_idx,
           frame_idx -
           (bptbl->oldest_bp == NO_BP
            ? 0 : garray_ent(bptbl->retired, bp_t, bptbl->oldest_bp).frame + 1)));
    sbevent_signal(bptbl->evt);
}

int
bptbl_push_frame(bptbl_t *bptbl, int oldest_bp)
{
    int frame_idx = bptbl->n_frame;

    E_DEBUG(2,("pushing frame %d, oldest bp %d in frame %d\n",
               frame_idx, oldest_bp,
               oldest_bp == NO_BP
               ? -1 : bptbl_ent(bptbl, oldest_bp)->frame));
    garray_expand_to(bptbl->ef_idx, frame_idx + 1);
    if (frame_idx - bptbl_active_frame(bptbl) >= bptbl->n_frame_alloc) {
        assert(bptbl->n_frame_alloc != 0);
        bptbl->n_frame_alloc *= 2;
        bptbl->valid_fr = bitvec_realloc(bptbl->valid_fr, bptbl->n_frame_alloc);
    }
    garray_ent(bptbl->ef_idx, bpidx_t, frame_idx) = bptbl_end_idx(bptbl);
    bptbl->n_frame = frame_idx + 1;
    bptbl_gc(bptbl, oldest_bp, frame_idx);
    return frame_idx;
}

int
bptbl_commit(bptbl_t *bptbl)
{
    bpidx_t src, dest, eidx;
    int frame_idx, dest_s_idx;

    /* This is the frame we're working on and its bps. */
    sbmtx_lock(bptbl->mtx);
    frame_idx = bptbl->n_frame - 1;
    dest = bptbl_ef_idx_internal(bptbl, frame_idx);
    dest_s_idx = garray_ent(bptbl->ent, bp_t, dest).s_idx;
    eidx = bptbl_end_idx(bptbl);
    E_DEBUG(4,("compacting %d bps\n", eidx - dest));
    if (eidx == dest) {
        sbmtx_unlock(bptbl->mtx);
        return 0;
    }
    /* Remove invalid bps and their rc entries. */
    for (src = dest; src < eidx; ++src) {
        bp_t *src_ent = garray_ptr(bptbl->ent, bp_t, src);
        bp_t *dest_ent = garray_ptr(bptbl->ent, bp_t, dest);
        if (src_ent->valid) {
            int rcsize = bptbl_rcsize(bptbl, src_ent);
            if (src_ent->s_idx != dest_s_idx) {
                E_DEBUG(4,("Moving %d rc scores from %d to %d for bptr %d\n",
                           rcsize, src_ent->s_idx, dest_s_idx, dest));
                assert(src_ent->s_idx > dest_s_idx);
                if (src < bptbl_end_idx(bptbl) - 1) {
                    bp_t *src1_ent = garray_ptr(bptbl->ent, bp_t, src + 1);
                    assert(dest_s_idx + rcsize <= src1_ent->s_idx);
                }
                /* Indices don't change but we still need to protect rc. */
                garray_move(bptbl->rc, dest_s_idx, src_ent->s_idx, rcsize);
            }
            *dest_ent = *src_ent;
            dest_ent->s_idx = dest_s_idx;
            dest_s_idx += rcsize;
            ++dest;
        }
    }
    E_DEBUG(4, ("Frame %d removed %d invalid bps out of %d\n",
                frame_idx, eidx - dest,
                eidx - bptbl_ef_idx_internal(bptbl, frame_idx)));

    /* Truncate active arrays. */
    garray_pop_from(bptbl->rc, dest_s_idx);
    garray_pop_from(bptbl->ent, dest);
    sbmtx_unlock(bptbl->mtx);

    return dest - src;
}

int
bptbl_finalize(bptbl_t *bptbl)
{
    int n_retired, first_retired_bp;

    sbmtx_lock(bptbl->mtx);
    E_DEBUG(2,("Final GC from frame %d to %d\n",
               bptbl_active_frame(bptbl), bptbl->n_frame));
    /* If there is nothing to GC then finish up. */
    if (bptbl_end_idx(bptbl) == bptbl_active_idx(bptbl)) {
        sbmtx_unlock(bptbl->mtx);
        return 0;
    }
    /* Expand the permutation table if necessary (probably). */
    garray_expand_to(bptbl->permute, bptbl_end_idx(bptbl));
    garray_set_base(bptbl->permute, bptbl_active_idx(bptbl));
    /* Mark and GC everything from the last frame. */
    n_retired = bptbl_mark(bptbl, bptbl->n_frame - 1, bptbl->n_frame);
    /* Include the last frame in the retired count. */
    n_retired += bptbl_ef_count_internal(bptbl, bptbl->n_frame - 1);
    E_DEBUG(2,("About to retire %d bps\n", n_retired));
    first_retired_bp = bptbl_retired_idx(bptbl);
    bptbl_retire(bptbl, n_retired, bptbl_end_idx(bptbl));
    bptbl_remap(bptbl, first_retired_bp,
                bptbl_end_idx(bptbl), bptbl_end_idx(bptbl));
    /* Just invalidate active entries, no need to move anything. */
    /* Empty the active entry table and set its base index to the
     * final index, which implicitly sets bptbl_active_idx (sorry
     * about that).  Do the same thing for bptbl_active_frame. */
    garray_reset_to(bptbl->ent, bptbl_end_idx(bptbl));
    garray_reset_to(bptbl->ef_idx, bptbl->n_frame);
    E_DEBUG(2,("Retired %d bps: now %d retired, %d active, first_active_sf %d\n",
               n_retired,
               bptbl_retired_idx(bptbl),
               bptbl_end_idx(bptbl) - bptbl_active_idx(bptbl),
               bptbl->oldest_bp == NO_BP
               ? 0 : garray_ent(bptbl->retired, bp_t, bptbl->oldest_bp).frame + 1));
    E_INFO("Allocated %d active and %d retired entries (%d + %d KiB)\n",
           garray_alloc_size(bptbl->ent),
           garray_alloc_size(bptbl->retired),
           garray_alloc_size(bptbl->ent) * sizeof(bp_t) / 1024,
           garray_alloc_size(bptbl->retired) * sizeof(bp_t) / 1024);
    E_INFO("Allocated %d right context scores (%d KiB)\n",
           garray_alloc_size(bptbl->rc),
           garray_alloc_size(bptbl->rc) * sizeof(int32) / 1024);
    E_INFO("Allocated %d permutation entries and %d end frame entries\n",
           garray_alloc_size(bptbl->permute), garray_alloc_size(bptbl->ef_idx));
    sbevent_signal(bptbl->evt);
    sbmtx_unlock(bptbl->mtx);
    return n_retired;
}

int
bptbl_release(bptbl_t *bptbl, bpidx_t first_idx)
{
    bpidx_t base_idx;
    bp_t *ent;

#if 0 /* For debugging purposes... */
    return 0;
#endif

    sbmtx_lock(bptbl->mtx);
    if (first_idx > bptbl_retired_idx(bptbl)) {
        E_DEBUG(2, ("%d outside retired, releasing up to %d\n",
                    first_idx, bptbl_retired_idx(bptbl)));
        first_idx = bptbl_retired_idx(bptbl);
    }

    base_idx = garray_base(bptbl->retired);
    E_DEBUG(2, ("Releasing bptrs from %d to %d\n",
                base_idx, first_idx));
    if (first_idx < base_idx) {
        sbmtx_unlock(bptbl->mtx);
        return 0;
    }

    ent = garray_ptr(bptbl->retired, bp_t, first_idx);
    garray_shift_from(bptbl->rc, ent->s_idx);
    garray_set_base(bptbl->rc, ent->s_idx);
    garray_shift_from(bptbl->retired, first_idx);
    garray_set_base(bptbl->retired, first_idx);

    sbmtx_unlock(bptbl->mtx);
    return first_idx - base_idx;
}

static bpidx_t
bptbl_find_exit(bptbl_t *bptbl, int32 wid)
{
    bpidx_t start, end, best;
    int32 best_score;
    int ef;

    if (bptbl_end_idx(bptbl) == 0)
        return NO_BP;

    /* We always take the last available frame, no matter what it
     * happens to be.  So take the last entry and scan backwards to
     * find the extents of its frame. */
    if (bptbl_active_idx(bptbl) == bptbl_end_idx(bptbl)) {
        bpidx_t first_retired = garray_base(bptbl->retired);
        /* Final, so it's in retired. */
        start = end = bptbl_retired_idx(bptbl);
        ef = bptbl_ent_internal(bptbl, start)->frame;
        while (start >= first_retired) {
            if (bptbl_ent_internal(bptbl, start)->frame != ef)
                break;
            --start;
        }
    }
    else {
        bpidx_t first_ent = garray_base(bptbl->ent);
        /* Not final, so it's in ent. */
        start = end = garray_next_idx(bptbl->ent) - 1;
        ef = bptbl_ent_internal(bptbl, start)->frame;
        while (start >= first_ent) {
            if (bptbl_ent_internal(bptbl, start)->frame != ef)
                break;
            --start;
        }
    }
    ++start;
    if (ef != bptbl->n_frame - 1) {
        E_WARN("No exits in final frame %d, using frame %d instead\n",
               bptbl->n_frame - 1, ef);
    }
    best_score = WORST_SCORE;
    best = NO_BP;
    while (start <= end) {
        bp_t *ent = bptbl_ent_internal(bptbl, start);
        if (ent->score BETTER_THAN best_score
            && (wid == BAD_S3WID || ent->wid == wid)) {
            best = start;
            best_score = ent->score;
        }
        ++start;
    }
    return best;
}

static int32
bptbl_ef_idx_internal(bptbl_t *bptbl, int frame_idx)
{
    if (frame_idx < bptbl_active_frame(bptbl))
        return 0;
    else if (frame_idx >= bptbl->n_frame)
        return bptbl_end_idx(bptbl);
    else
        return garray_ent(bptbl->ef_idx, bpidx_t, frame_idx);
}

int32
bptbl_ef_idx(bptbl_t *bptbl, int frame_idx)
{
    int32 idx;

    sbmtx_lock(bptbl->mtx);
    idx = bptbl_ef_idx_internal(bptbl, frame_idx);
    sbmtx_unlock(bptbl->mtx);

    return idx;
}

static bp_t *
bptbl_ent_internal(bptbl_t *bptbl, bpidx_t bpidx)
{
    if (bpidx == NO_BP)
        return NULL;
    if (bpidx < bptbl_active_idx(bptbl))
        return garray_ptr(bptbl->retired, bp_t, bpidx);
    else
        return garray_ptr(bptbl->ent, bp_t, bpidx);
}

int
bptbl_get_bp(bptbl_t *bptbl, bpidx_t bpidx, bp_t *out_bp)
{
    bp_t *ent;

    if (bpidx == NO_BP)
        return -1;

    sbmtx_lock(bptbl->mtx);
    ent = bptbl_ent_internal(bptbl, bpidx);
    if (ent == NULL) {
        sbmtx_unlock(bptbl->mtx);
        return -1;
    }
    memcpy(out_bp, ent, sizeof(*out_bp));
    sbmtx_unlock(bptbl->mtx);

    return 0;
}

int
bptbl_set_bp(bptbl_t *bptbl, bpidx_t bpidx, bp_t const *bp)
{
    bp_t *ent;

    if (bpidx == NO_BP)
        return -1;

    sbmtx_lock(bptbl->mtx);
    ent = bptbl_ent_internal(bptbl, bpidx);
    if (ent == NULL) {
        sbmtx_unlock(bptbl->mtx);
        return -1;
    }
    memcpy(ent, bp, sizeof(*ent));
    sbmtx_unlock(bptbl->mtx);

    return 0;
}

bpidx_t
bptbl_active_idx(bptbl_t *bptbl)
{
    return garray_base(bptbl->ent);
}

bpidx_t
bptbl_retired_idx(bptbl_t *bptbl)
{
    return garray_next_idx(bptbl->retired);
}


bpidx_t
bptbl_end_idx(bptbl_t *bptbl)
{
    return garray_next_idx(bptbl->ent);
}

int
bptbl_active_frame(bptbl_t *bptbl)
{
    return garray_base(bptbl->ef_idx);
}

int
bptbl_frame_idx(bptbl_t *bptbl)
{
    return bptbl->n_frame;
}

int
bptbl_active_sf(bptbl_t *bptbl)
{
    int sf;
    sbmtx_lock(bptbl->mtx);
    if (bptbl->oldest_bp == NO_BP)
        sf = 0;
    else
        sf = garray_ent(bptbl->ent, bp_t, bptbl->oldest_bp).frame + 1;
    sbmtx_unlock(bptbl->mtx);
    return sf;
}

int
bptbl_wait(bptbl_t *bptbl, int timeout)
{
    int s = (timeout == -1) ? -1 : 0;
    int rc;

    if ((rc = sbevent_wait(bptbl->evt, s, timeout)) < 0)
        return rc;
    return 0;
}

static int
bptbl_sf_internal(bptbl_t *bptbl, bpidx_t bpidx)
{
    bp_t *ent;
    bp_t *prev;

    ent = bptbl_ent_internal(bptbl, bpidx);
    if (ent == NULL)
        return -1;
    else {
        prev = bptbl_ent_internal(bptbl, ent->bp);
        if (prev == NULL)
            return 0;
        else
            return prev->frame + 1;
    }
}

int
bptbl_sf(bptbl_t *bptbl, bpidx_t bpidx)
{
    int sf;

    sbmtx_lock(bptbl->mtx);
    sf = bptbl_sf_internal(bptbl, bpidx);
    sbmtx_unlock(bptbl->mtx);

    return sf;
}

static int
bptbl_ef_count_internal(bptbl_t *bptbl, int frame_idx)
{
    bpidx_t start, end;
    start = bptbl_ef_idx_internal(bptbl, frame_idx);
    end = bptbl_ef_idx_internal(bptbl, frame_idx + 1);
    return end - start;
}

int
bptbl_ef_count(bptbl_t *bptbl, int frame_idx)
{
    int count;

    sbmtx_lock(bptbl->mtx);
    count = bptbl_ef_count_internal(bptbl, frame_idx);
    sbmtx_unlock(bptbl->mtx);

    return count;
}

void
bptbl_set_rcscore(bptbl_t *bptbl, bpidx_t bpidx, int rc, int32 score)
{
    bp_t *bpe;

    sbmtx_lock(bptbl->mtx);
    bpe = bptbl_ent_internal(bptbl, bpidx);
    garray_ent(bptbl->rc, int32, bpe->s_idx + rc) = score;
    sbmtx_unlock(bptbl->mtx);
}

int
bptbl_get_rcscores(bptbl_t *bptbl, bpidx_t bpidx, int32 *out_rcscores)
{
    bp_t *bpe;
    int rcsize;

    sbmtx_lock(bptbl->mtx);
    bpe = bptbl_ent_internal(bptbl, bpidx);
    rcsize = bptbl_rcsize(bptbl, bpe);
    memcpy(out_rcscores,
           garray_ptr(bptbl->rc, int32, bpe->s_idx),
           rcsize * sizeof(int32));
    sbmtx_unlock(bptbl->mtx);
    return rcsize;
}

/**
 * Get the number of right context entries for a backpointer.
 *
 * This is only used internally, so it has no locking.
 */
static int
bptbl_rcsize(bptbl_t *bptbl, bp_t *be)
{
    int rcsize;

    if (dict_is_single_phone(bptbl->d2p->dict, be->wid)) {
        be->last2_phone = -1;
        rcsize = 1;
    }
    else {
        be->last2_phone = dict_second_last_phone(bptbl->d2p->dict, be->wid);
        rcsize = dict2pid_rssid(bptbl->d2p, be->last_phone, be->last2_phone)->n_ssid;
    }

    return rcsize;
}

static void
bptbl_fake_lmstate_internal(bptbl_t *bptbl, bp_t *ent)
{
    bp_t *prev;

    prev = bptbl_ent_internal(bptbl, ent->bp);
    /* Propagate lm state for fillers, rotate it for words. */
    if (dict_filler_word(bptbl->d2p->dict, ent->wid)) {
        if (prev != NULL) {
            ent->real_wid = prev->real_wid;
            ent->prev_real_wid = prev->prev_real_wid;
        }
    }
    else {
        ent->real_wid = dict_basewid(bptbl->d2p->dict, ent->wid);
        if (prev != NULL)
            ent->prev_real_wid = prev->real_wid;
    }
}

bpidx_t
bptbl_enter(bptbl_t *bptbl, int32 w, int32 path, int32 score, int rc)
{
    int32 i, rcsize, *bss;
    size_t bss_head;
    bpidx_t bpidx;
    bp_t be, *bpe;

    /* This might happen if recognition fails. */
    if (bptbl_end_idx(bptbl) == NO_BP) {
        E_ERROR("No entries in backpointer table!");
        return NO_BP;
    }

    sbmtx_lock(bptbl->mtx);
    /* Append a new backpointer. */
    memset(&be, 0, sizeof(be));
    be.wid = w;
    be.frame = bptbl->n_frame - 1;
    be.bp = path;
    be.score = score;
    be.s_idx = garray_next_idx(bptbl->rc);
    be.valid = TRUE;
    be.last_phone = dict_last_phone(bptbl->d2p->dict, w);
    bpidx = garray_next_idx(bptbl->ent);
    bpe = garray_append(bptbl->ent, &be);
    /* Set up its LM state */
    bptbl_fake_lmstate_internal(bptbl, bpe);

    /* DICT2PID */
    /* Get diphone ID for final phone and number of ssids corresponding to it. */
    rcsize = bptbl_rcsize(bptbl, bpe);
    /* Allocate some space on the bptbl->bscore_stack for all of these triphones. */
    /* Expand the bss table if necessary. */
    bss_head = be.s_idx;
    garray_expand_to(bptbl->rc, bss_head + rcsize);
    bss = garray_ptr(bptbl->rc, int32, bss_head);
      for (i = 0; i < rcsize; ++i)
        *bss++ = WORST_SCORE;
    garray_ent(bptbl->rc, int32, bss_head + rc) = score;

    E_DEBUG(3,("Entered bp %d sf %d ef %d s_idx %d active_fr %d\n",
               bptbl_end_idx(bptbl) - 1,
               bptbl_sf_internal(bptbl, bptbl_end_idx(bptbl) - 1),
               bptbl->n_frame - 1, bss_head, bptbl_active_frame(bptbl)));
    assert(bptbl_sf_internal(bptbl, bptbl_end_idx(bptbl) - 1)
           >= bptbl_active_frame(bptbl));
    sbmtx_unlock(bptbl->mtx);
    return bpidx;
}

void
bptbl_fake_lmstate(bptbl_t *bptbl, int32 bp, int32 new_prev)
{
    bp_t *ent;

    assert(bp != NO_BP);
    sbmtx_lock(bptbl->mtx);
    ent = bptbl_ent_internal(bptbl, bp);
    ent->bp = new_prev;
    bptbl_fake_lmstate_internal(bptbl, ent);
    sbmtx_unlock(bptbl->mtx);
}

char *
bptbl_hyp(bptbl_t *bptbl, int32 *out_score, int32 finish_wid)
{
    bpidx_t exit;
    bp_t *bpe;
    char *c, *hyp_str;
    size_t len;

    sbmtx_lock(bptbl->mtx);
    /* Look for </s> in the last frame. */
    if ((exit = bptbl_find_exit(bptbl, finish_wid)) == NO_BP) {
        /* If not found then take the best scoring word exit (with warning). */
        exit = bptbl_find_exit(bptbl, BAD_S3WID);
        if (exit == NO_BP) {
            E_ERROR("No word exits in last frame: recognition failure?\n");
            sbmtx_unlock(bptbl->mtx);
            return NULL;
        }
        E_WARN("No %s found in last frame, using %s instead\n",
               dict_wordstr(bptbl->d2p->dict, finish_wid),
               dict_wordstr(bptbl->d2p->dict,
                            bptbl_ent_internal(bptbl, exit)->wid));
    }

    bpe = bptbl_ent_internal(bptbl, exit);
    if (out_score)
        *out_score = bpe->score;
    len = 0;
    while (bpe != NULL) {
        assert(bpe->valid);
        if (dict_real_word(bptbl->d2p->dict, bpe->wid))
            len += strlen(dict_basestr(bptbl->d2p->dict, bpe->wid)) + 1;
        bpe = bptbl_ent_internal(bptbl, bpe->bp);
    }

    if (len == 0) {
        sbmtx_unlock(bptbl->mtx);
	return NULL;
    }
    hyp_str = ckd_calloc(1, len);

    bpe = bptbl_ent_internal(bptbl, exit);
    c = hyp_str + len - 1;
    while (bpe != NULL) {
        size_t len;
        if (dict_real_word(bptbl->d2p->dict, bpe->wid)) {
            len = strlen(dict_basestr(bptbl->d2p->dict, bpe->wid));
            c -= len;
            memcpy(c, dict_basestr(bptbl->d2p->dict, bpe->wid), len);
            if (c > hyp_str) {
                --c;
                *c = ' ';
            }
        }
        bpe = bptbl_ent_internal(bptbl, bpe->bp);
    }

    sbmtx_unlock(bptbl->mtx);
    return hyp_str;
}

static void
bptbl_bp2itor(ps_seg_t *seg, int bp)
{
    bptbl_seg_t *bseg = (bptbl_seg_t *)seg;
    bp_t *be, *pbe;

    sbmtx_lock(bseg->bptbl->mtx);
    be = bptbl_ent_internal(bseg->bptbl, bp);
    pbe = bptbl_ent_internal(bseg->bptbl, be->bp);
    seg->word = dict_wordstr(bseg->bptbl->d2p->dict, be->wid);
    seg->ef = be->frame;
    seg->sf = pbe ? pbe->frame + 1 : 0;
    seg->prob = 0; /* Bogus value... */
    /* Compute acoustic and LM scores for this segment. */
    if (pbe == NULL) {
        seg->ascr = be->score;
        seg->lscr = 0;
        seg->lback = 0;
    }
    else {
        /* Language model score calculation (for what it's worth,
         * which isn't much since we don't have real language model
         * state) is search-dependent. */
        seg->ascr = be->score - pbe->score;
        seg->lscr = 0;
    }
    sbmtx_unlock(bseg->bptbl->mtx);
}

static void
bptbl_seg_free(ps_seg_t *seg)
{
    bptbl_seg_t *itor = (bptbl_seg_t *)seg;
    
    ckd_free(itor->bpidx);
    ckd_free(itor);
}

static ps_seg_t *
bptbl_seg_next(ps_seg_t *seg)
{
    bptbl_seg_t *itor = (bptbl_seg_t *)seg;

    if (++itor->cur == itor->n_bpidx) {
        bptbl_seg_free(seg);
        return NULL;
    }

    bptbl_bp2itor(seg, itor->bpidx[itor->cur]);
    return seg;
}

static ps_segfuncs_t bptbl_segfuncs = {
    /* seg_next */ bptbl_seg_next,
    /* seg_free */ bptbl_seg_free
};

ps_seg_t *
bptbl_seg_iter(bptbl_t *bptbl, int32 *out_score, int32 finish_wid)
{
    bptbl_seg_t *itor;
    bpidx_t exit;
    bp_t *bpe;
    int cur;

    sbmtx_lock(bptbl->mtx);
    /* Look for </s> in the last frame. */
    if ((exit = bptbl_find_exit(bptbl, finish_wid)) == NO_BP) {
        /* If not found then take the best scoring word exit (with warning). */
        exit = bptbl_find_exit(bptbl, BAD_S3WID);
        if (exit == NO_BP) {
            E_ERROR("No word exits in last frame: recognition failure?\n");
            sbmtx_unlock(bptbl->mtx);
            return NULL;
        }
        E_WARN("No %s found in last frame, using %s instead\n",
               dict_wordstr(bptbl->d2p->dict, finish_wid),
               dict_wordstr(bptbl->d2p->dict,
                            bptbl_ent_internal(bptbl, exit)->wid));
    }

    bpe = bptbl_ent_internal(bptbl, exit);
    if (out_score)
        *out_score = bpe->score;

    /* Calling this an "iterator" is a bit of a misnomer since we have
     * to get the entire backtrace in order to produce it.  On the
     * other hand, all we actually need is the bptbl IDs, and we can
     * allocate a fixed-size array of them. */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->base.vt = &bptbl_segfuncs;
    itor->base.lwf = 1.0;
    itor->bptbl = bptbl;
    itor->n_bpidx = 0;
    while (bpe != NULL) {
        ++itor->n_bpidx;
        bpe = bptbl_ent_internal(bptbl, bpe->bp);
    }
    if (itor->n_bpidx == 0) {
        ckd_free(itor);
        sbmtx_unlock(bptbl->mtx);
        return NULL;
    }
    itor->bpidx = ckd_calloc(itor->n_bpidx, sizeof(*itor->bpidx));
    cur = itor->n_bpidx - 1;
    bpe = bptbl_ent_internal(bptbl, exit);
    while (bpe != NULL) {
        itor->bpidx[cur] = exit;
        exit = bpe->bp;
        bpe = bptbl_ent_internal(bptbl, exit);
        --cur;
    }

    sbmtx_unlock(bptbl->mtx);

    /* Fill in relevant fields for first element. */
    bptbl_bp2itor((ps_seg_t *)itor, itor->bpidx[0]);

    return (ps_seg_t *)itor;
}
