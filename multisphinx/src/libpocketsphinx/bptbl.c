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

static int bptbl_rcsize(bptbl_t *bptbl, bp_t *be);

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

    bptbl->n_permute_alloc = n_frame_alloc;
    bptbl->permute = ckd_calloc(bptbl->n_permute_alloc, sizeof(*bptbl->permute));
    E_INFO("Allocated %d KiB for permutation table\n",
           bptbl->n_permute_alloc * sizeof(*bptbl->permute) / 1024);

    bptbl->bscore_stack_size = n_alloc * 20;
    bptbl->bscore_stack = ckd_calloc(bptbl->bscore_stack_size,
                                     sizeof(*bptbl->bscore_stack));
    E_INFO("Allocated %d KiB for right context scores\n",
           bptbl->bscore_stack_size * sizeof(*bptbl->bscore_stack) / 1024);

    bptbl->n_frame_alloc = n_frame_alloc;
    bptbl->ef_idx = ckd_calloc(bptbl->n_frame_alloc,
                               sizeof(*bptbl->ef_idx));
    bptbl->valid_fr = bitvec_alloc(bptbl->n_frame_alloc);

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
    ckd_free(bptbl->permute);
    ckd_free(bptbl->bscore_stack);
    ckd_free(bptbl->ef_idx);
    bitvec_free(bptbl->valid_fr);
    ckd_free(bptbl);
    return 0;
}

void
bptbl_reset(bptbl_t *bptbl)
{
    int i;

    for (i = 0; i < bptbl->n_frame_alloc; ++i) {
        bptbl->ef_idx[i] = -1;
    }
    bitvec_clear_all(bptbl->valid_fr, bptbl->n_frame_alloc);
    garray_reset(bptbl->ent);
    garray_reset(bptbl->retired);
    bptbl->first_invert_bp = 0;
    bptbl->dest_s_idx = 0;
    bptbl->n_frame = 0;
    bptbl->bss_head = 0;
    bptbl->active_fr = 0;
    bptbl->oldest_bp = NO_BP;
}

void
bptbl_dump(bptbl_t *bptbl)
{
    int i;

    E_INFO("Retired backpointers (%d entries, oldest active %d):\n",
           bptbl_retired_idx(bptbl), bptbl->oldest_bp);
    for (i = 0; i < bptbl_retired_idx(bptbl); ++i) {
        bp_t *ent = bptbl_ent(bptbl, i);
        assert(ent->valid);
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                    i, dict_wordstr(bptbl->d2p->dict, ent->wid),
                    bptbl_sf(bptbl, i),
                    ent->frame,
                    ent->score,
                    ent->bp);
    }
    E_INFO("Active backpointers (%d entries starting at %d):\n",
           bptbl_end_idx(bptbl) - bptbl->ef_idx[0], bptbl->ef_idx[0]);
    for (i = bptbl->ef_idx[0]; i < bptbl_end_idx(bptbl); ++i) {
        bp_t *ent = bptbl_ent(bptbl, i);
        if (!ent->valid)
            E_INFO_NOFN("%-5d INVALID\n", i);
        else
            E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                        i, dict_wordstr(bptbl->d2p->dict, ent->wid),
                        bptbl_sf(bptbl, i),
                        ent->frame,
                        ent->score,
                        ent->bp);
    }
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
    assert(ef > bptbl->active_fr);

    /* Invalidate all active backpointer entries up to ef. */
    E_DEBUG(2,("Invalidating backpointers from %d to %d (%d to %d)\n",
               bptbl->active_fr, ef,
               bptbl_ef_idx(bptbl, bptbl->active_fr), bptbl_ef_idx(bptbl, ef)));
    for (i = bptbl_ef_idx(bptbl, bptbl->active_fr);
         i < bptbl_ef_idx(bptbl, ef); ++i) {
        E_DEBUG(5,("Invalidate bp %d\n", i));
        garray_ent(bptbl->ent, bp_t, i).valid = FALSE;
    }

    /* Now re-activate all ones backwards reachable from the search graph. */
    E_DEBUG(2,("Finding coaccessible frames from backpointers from %d to %d (%d to %d)\n",
               ef, cf,
               bptbl_ef_idx(bptbl, ef), bptbl_ef_idx(bptbl, cf)));
    /* Mark everything immediately reachable from (ef..cf) */
    bitvec_clear_all(bptbl->valid_fr, cf - bptbl->active_fr);
    n_active_fr = 0;
    /* NOTE: This for statement can be sped up at the cost of being
     * less obvious. */
    for (i = bptbl_ef_idx(bptbl, ef);
         i < bptbl_ef_idx(bptbl, cf); ++i) {
        bp_t *ent, *prev;
        ent = bptbl_ent(bptbl, i);
        prev = bptbl_prev(bptbl, ent);
        assert(ent->valid);
        if (prev != NULL) {
            int frame = prev->frame;
            if (frame >= bptbl->active_fr
                && bitvec_is_clear(bptbl->valid_fr, frame - bptbl->active_fr)) {
                E_DEBUG(5,("Validate frame %d\n", frame - bptbl->active_fr));
                bitvec_set(bptbl->valid_fr, frame - bptbl->active_fr);
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
        for (i = bptbl->active_fr; i <= last_gc_fr; ++i) {
            /* NOTE: We asserted i >= bptbl->active_fr */
            if (bitvec_is_set(bptbl->valid_fr, i - bptbl->active_fr)) {
                bitvec_clear(bptbl->valid_fr, i - bptbl->active_fr);
                /* Add all backpointers in this frame (the bogus
                 * lattice generation algorithm) */
                for (j = bptbl_ef_idx(bptbl, i);
                     j < bptbl_ef_idx(bptbl, i + 1); ++j) {
                    bp_t *ent = bptbl_ent(bptbl, j);
                    bp_t *prev = bptbl_prev(bptbl, ent);
                    ent->valid = TRUE;
                    if (prev != NULL) {
                        int frame = prev->frame;
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
    for (j = 0, i = bptbl_ef_idx(bptbl, bptbl->active_fr);
         i < bptbl_ef_idx(bptbl, ef); ++i) {
        bp_t *ent = bptbl_ent(bptbl, i);
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
    dest = bptbl->first_invert_bp;
    /* Expand retired if necessary. */
    garray_expand(bptbl->retired, dest + n_retired);
    /* Note we use the "cooked" backpointer indices here. */
    for (src = bptbl->ef_idx[0]; src < eidx; ++src) {
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
                memmove(bptbl->bscore_stack + bptbl->dest_s_idx,
                        bptbl->bscore_stack + src_ent->s_idx,
                        rcsize * sizeof(*bptbl->bscore_stack));
            }
            *dest_ent = *src_ent;
            dest_ent->s_idx = bptbl->dest_s_idx;

            /* FIXME: Going to move permute over to garray + base_idx too... */
            assert(src - bptbl->ef_idx[0] < bptbl->n_permute_alloc);
            bptbl->permute[src - bptbl->ef_idx[0]] = dest;
            bptbl->dest_s_idx += rcsize;
            ++dest;
        }
        else {
            E_DEBUG(4,("permute %d => -1 src %d ef_idx %d\n",
                       src, src - bptbl->ef_idx[0], bptbl->ef_idx[0]));
            assert(src - bptbl->ef_idx[0] < bptbl->n_permute_alloc);
            bptbl->permute[src - bptbl->ef_idx[0]] = -1;
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
            memmove(bptbl->bscore_stack + active_dest_s_idx,
                    bptbl->bscore_stack + src_ent->s_idx,
                    rcsize * sizeof(*bptbl->bscore_stack));
            src_ent->s_idx = active_dest_s_idx;
            active_dest_s_idx += rcsize;
            ++src;
        }
        bptbl->bss_head = active_dest_s_idx;
    }
    return dest;
}


/**
 * Remap backpointer indices in backpointer table.
 *
 * This needs to be done both in newly retired backpointers (which
 * range from first_invert_bp to last_retired_bp) as well as in still
 * active backpointers, which range from first_active_bp to bptbl->n_ent.
 *
 * This has the side effect of setting bptbl->oldest_bp
 *
 * @param bptbl Backpointer table.
 * @param last_retired_bp Index of last retired backpointer plus one
 * @param last_remapped_bp Last backpointer index to be remapped
 * @param first_active_bp Index of first remaining active backpointer
 */
static void
bptbl_remap(bptbl_t *bptbl, int last_retired_bp,
            int last_remapped_bp, int first_active_bp)
{
    int i;

    E_DEBUG(2,("inverting %d:%d from %d to %d and %d to %d\n",
               bptbl->first_invert_bp, last_remapped_bp,
               bptbl->first_invert_bp, last_retired_bp,
               first_active_bp, bptbl_end_idx(bptbl)));
    /* First remap backpointers in newly retired bps. */
    for (i = bptbl->first_invert_bp; i < last_retired_bp; ++i) {
        bp_t *bpe = garray_ptr(bptbl->retired, bp_t, i);
        /* Remember, these are the *source* backpointer indices, so
         * they fall in the range between prev_active_fr (which is the
         * first index of ef_idx) and active_fr. */
        if (bpe->bp >= bptbl->ef_idx[0]
            && bpe->bp < last_remapped_bp) {
            assert(bpe->bp - bptbl->ef_idx[0] < bptbl->n_permute_alloc);
            if (bpe->bp != bptbl->permute[bpe->bp - bptbl->ef_idx[0]])
                E_DEBUG(4,("remap retired %d => %d in %d\n",
                           bpe->bp, bptbl->permute[bpe->bp - bptbl->ef_idx[0]], i));
            bpe->bp = bptbl->permute[bpe->bp - bptbl->ef_idx[0]];
            assert(bptbl_sf(bptbl, i) <= bpe->frame);
        }
    }
    /* Now remap backpointers in still-active bps (which point to the
     * newly retired ones) */
    bptbl->oldest_bp = last_retired_bp - 1;
    for (i = first_active_bp; i < bptbl_end_idx(bptbl); ++i) {
        bp_t *bpe = garray_ptr(bptbl->ent, bp_t, i);
        if (bpe->bp >= bptbl->ef_idx[0] && bpe->bp < last_remapped_bp) {
            assert(bpe->bp - bptbl->ef_idx[0] < bptbl->n_permute_alloc);
            if (bpe->bp != bptbl->permute[bpe->bp - bptbl->ef_idx[0]])
                E_DEBUG(4,("remap active %d => %d in %d\n",
                           bpe->bp, bptbl->permute[bpe->bp - bptbl->ef_idx[0]], i));
            bpe->bp = bptbl->permute[bpe->bp - bptbl->ef_idx[0]];
            assert(bptbl_sf(bptbl, i) <= bpe->frame);
        }
        if (bpe->bp < bptbl->oldest_bp)
            bptbl->oldest_bp = bpe->bp;
    }
}

/**
 * Update the active frame pointer and backpointer array.
 */
static void
bptbl_update_active(bptbl_t *bptbl, int active_fr, int last_retired_bp)
{
    int frame_delta = active_fr - bptbl->active_fr;
    int bp_delta;
    /* This means nothing happened. */
    if (frame_delta == 0)
        return;
    assert(frame_delta > 0);
    bp_delta = bptbl->ef_idx[frame_delta] - bptbl->ef_idx[0];

    /* Push back active backpointers (eventually this will be circular) */
    E_DEBUG(3,("moving %d ent from %d (%d - %d)\n",
               bptbl_end_idx(bptbl) - bptbl->ef_idx[frame_delta],
               bp_delta, bptbl->ef_idx[frame_delta], bptbl->ef_idx[0]));
    garray_shift(bptbl->ent, bp_delta);
    
    /* Update ef_idx (implicitly updating ef_idx[0] */
    E_DEBUG(3,("moving %d ef_idx from %d (%d - %d)\n",
               bptbl->n_frame - active_fr,
               frame_delta, active_fr, bptbl->active_fr));
    assert(frame_delta + bptbl->n_frame - active_fr < bptbl->n_frame_alloc);
    memmove(bptbl->ef_idx, bptbl->ef_idx + frame_delta,
            (bptbl->n_frame - active_fr) * sizeof(*bptbl->ef_idx));

    /* And now update stuff. */
    garray_set_base(bptbl->ent, bptbl->ef_idx[0]);
    bptbl->active_fr = active_fr;
    bptbl->first_invert_bp = last_retired_bp;
}

static void
bptbl_gc(bptbl_t *bptbl, int oldest_bp, int frame_idx)
{
    int prev_active_fr, active_fr;
    int n_retired, last_retired_bp;

    /* active_fr is the first frame which is still active in search
     * (i.e. for which outgoing word arcs can still be generated).
     * Therefore, any future backpointer table entries will not point
     * backwards to any backpointers before active_fr, and thus any
     * backpointers which are not reachable from those exiting in
     * active_fr will never be reachable. */
    prev_active_fr = bptbl->active_fr;
    if (oldest_bp == NO_BP)
        active_fr = 0;
    else
        active_fr = bptbl_ent(bptbl, oldest_bp)->frame;
    assert(active_fr >= prev_active_fr);
    /* Need at least 2 frames to GC. */
    if (active_fr <= prev_active_fr + 1)
        return;
    /* If there is nothing to GC then finish up. */
    if (bptbl_ef_idx(bptbl, prev_active_fr)
        == bptbl_ef_idx(bptbl, active_fr)) {
        bptbl_update_active(bptbl, active_fr, bptbl->first_invert_bp);
        return;
    }
    E_DEBUG(2,("GC from frame %d to %d\n", prev_active_fr, active_fr));
    /* Expand the permutation table if necessary. */
    if ((bptbl_ef_idx(bptbl, active_fr) - bptbl->ef_idx[0])
        > bptbl->n_permute_alloc) {
        while (bptbl_ef_idx(bptbl, active_fr) - bptbl->ef_idx[0]
               > bptbl->n_permute_alloc)
            bptbl->n_permute_alloc *= 2;
        bptbl->permute = ckd_realloc(bptbl->permute,
                                     bptbl->n_permute_alloc
                                     * sizeof(*bptbl->permute));
        E_INFO("Resized permutation table to %d entries (%d KiB)\n",
               bptbl->n_permute_alloc,
               bptbl->n_permute_alloc * sizeof(*bptbl->permute) / 1024);
    }
    /* Mark, compact, snap pointers. */
    n_retired = bptbl_mark(bptbl, active_fr, frame_idx);
    E_DEBUG(2,("About to retire %d bps\n", n_retired));
    last_retired_bp = bptbl_retire(bptbl, n_retired, bptbl_ef_idx(bptbl, active_fr));
    assert(n_retired == last_retired_bp - bptbl->first_invert_bp);
    bptbl_remap(bptbl, last_retired_bp,
                bptbl_ef_idx(bptbl, active_fr),
                bptbl_ef_idx(bptbl, active_fr));
    bptbl_update_active(bptbl, active_fr, last_retired_bp);
    E_INFO("Retired %d bps: now %d retired, %d active\n", n_retired,
           bptbl_retired_idx(bptbl),
           bptbl_end_idx(bptbl) - bptbl->ef_idx[0]);
    E_INFO("First active sf %d output frame %d window %d\n",
           bptbl->oldest_bp == NO_BP
           ? 0 : garray_ent(bptbl->retired, bp_t, bptbl->oldest_bp).frame + 1,
           frame_idx,
           frame_idx -
           (bptbl->oldest_bp == NO_BP
            ? 0 : garray_ent(bptbl->retired, bp_t, bptbl->oldest_bp).frame + 1));
}

int
bptbl_push_frame(bptbl_t *bptbl, int oldest_bp)
{
    int frame_idx = bptbl->n_frame;

    E_DEBUG(2,("pushing frame %d, oldest bp %d in frame %d\n",
               frame_idx, oldest_bp,
               oldest_bp == NO_BP
               ? -1 : bptbl_ent(bptbl, oldest_bp)->frame));
    if (frame_idx - bptbl->active_fr >= bptbl->n_frame_alloc) {
        bptbl->n_frame_alloc *= 2;
        E_INFO("Reallocating frame-based bptr arrays to %d\n", bptbl->n_frame_alloc);
        bptbl->ef_idx = ckd_realloc(bptbl->ef_idx,
                                    bptbl->n_frame_alloc * sizeof(*bptbl->ef_idx));
        bptbl->valid_fr = bitvec_realloc(bptbl->valid_fr, bptbl->n_frame_alloc);
    }
    bptbl->ef_idx[frame_idx - bptbl->active_fr] = bptbl_end_idx(bptbl);
    bptbl->n_frame = frame_idx + 1;
    bptbl_gc(bptbl, oldest_bp, frame_idx);
    return frame_idx;
}

int
bptbl_finalize(bptbl_t *bptbl)
{
    int n_retired, last_retired_bp;

    E_DEBUG(2,("Final GC from frame %d to %d\n",
               bptbl->active_fr, bptbl->n_frame));
    /* If there is nothing to GC then finish up. */
    if (bptbl_end_idx(bptbl) == bptbl->ef_idx[0])
        return 0;
    /* Expand the permutation table if necessary (probably). */
    if (bptbl_end_idx(bptbl) - bptbl->ef_idx[0]
        > bptbl->n_permute_alloc) {
        while (bptbl_end_idx(bptbl) - bptbl->ef_idx[0]
               > bptbl->n_permute_alloc)
            bptbl->n_permute_alloc *= 2;
        bptbl->permute = ckd_realloc(bptbl->permute,
                                     bptbl->n_permute_alloc
                                     * sizeof(*bptbl->permute));
        E_INFO("Resized permutation table to %d entries (%d KiB)\n",
               bptbl->n_permute_alloc,
               bptbl->n_permute_alloc * sizeof(*bptbl->permute) / 1024);
    }
    /* Mark and GC everything from the last frame. */
    n_retired = bptbl_mark(bptbl, bptbl->n_frame - 1, bptbl->n_frame);
    /* Include the last frame in the retired count. */
    n_retired += bptbl_ef_count(bptbl, bptbl->n_frame - 1);
    E_DEBUG(2,("About to retire %d bps\n", n_retired));
    last_retired_bp = bptbl_retire(bptbl, n_retired, bptbl_end_idx(bptbl));
    bptbl_remap(bptbl, last_retired_bp, bptbl_end_idx(bptbl), bptbl_end_idx(bptbl));
    /* Just invalidate active entries, no need to move anything. */
    bptbl->first_invert_bp = last_retired_bp;
    bptbl->active_fr = bptbl->n_frame;
    /* FIXME: Will have to change this eventually when we use garray_t for ef_idx */
    bptbl->ef_idx[0] = bptbl_end_idx(bptbl);
    E_INFO("Retired %d bps: now %d retired, %d active, first_active_sf %d\n", n_retired,
           bptbl_retired_idx(bptbl),
           bptbl_end_idx(bptbl) - bptbl->ef_idx[0],
           bptbl->oldest_bp == NO_BP
           ? 0 : garray_ent(bptbl->retired, bp_t, bptbl->oldest_bp).frame + 1);
    E_INFO("Allocated %d active and %d retired entries\n",
           garray_alloc_size(bptbl->ent),
           garray_alloc_size(bptbl->retired));
    return n_retired;
}

bp_t *
bptbl_find_exit(bptbl_t *bptbl, int32 wid)
{
    bp_t *start, *end, *best;
    int32 best_score;
    int ef;

    if (bptbl_end_idx(bptbl) == 0)
        return NULL;

    /* We always take the last available frame, no matter what it
     * happens to be.  So take the last entry and scan backwards to
     * find the extents of its frame. */
    if (bptbl->ef_idx[0] == bptbl_end_idx(bptbl)) {
        bp_t *first_retired = garray_ptr(bptbl->retired, bp_t, 0);
        /* Final, so it's in retired. */
        start = garray_ptr(bptbl->retired, bp_t, bptbl->first_invert_bp - 1);
        end = garray_ptr(bptbl->retired, bp_t, bptbl->first_invert_bp);
        ef = start->frame;
        while (start >= first_retired && start->frame == ef)
            --start;
    }
    else {
        bp_t *first_ent = garray_ptr(bptbl->ent, bp_t, garray_base(bptbl->ent));
        /* Not final, so it's in ent. */
        start = garray_ptr(bptbl->ent, bp_t, garray_next_idx(bptbl->ent) - 1);
        end = garray_ptr(bptbl->ent, bp_t, garray_next_idx(bptbl->ent));
        ef = start->frame;
        /* FIXME: When bptbl->ent is circular this will no longer work. */
        while (start >= first_ent && start->frame == ef)
            --start;
    }
    ++start;
    if (ef != bptbl->n_frame - 1) {
        E_WARN("No exits in final frame %d, using frame %d instead\n",
               bptbl->n_frame - 1, ef);
    }
    E_INFO("Last frame %d has %d entries\n", ef, end - start);
    best_score = start->score;
    best = NULL;
    while (start < end) {
        if (start->score BETTER_THAN best_score
            && (wid == BAD_S3WID || start->wid == wid)) {
            best = start;
            best_score = best->score;
        }
        ++start;
    }
    return best;
}

int32
bptbl_ef_idx(bptbl_t *bptbl, int frame_idx)
{
    if (frame_idx < bptbl->active_fr)
        return 0;
    else if (frame_idx >= bptbl->n_frame)
        return bptbl_end_idx(bptbl);
    else {
        return bptbl->ef_idx[frame_idx - bptbl->active_fr];
    }
}

bp_t *
bptbl_ent(bptbl_t *bptbl, bpidx_t bpidx)
{
    if (bpidx == NO_BP)
        return NULL;
    if (bpidx < bptbl->ef_idx[0])
        return garray_ptr(bptbl->retired, bp_t, bpidx);
    else
        return garray_ptr(bptbl->ent, bp_t, bpidx);
}

bpidx_t
bptbl_active_idx(bptbl_t *bptbl)
{
    return garray_base(bptbl->ent);
}

bpidx_t
bptbl_retired_idx(bptbl_t *bptbl)
{
    return bptbl->first_invert_bp;
}


bpidx_t
bptbl_end_idx(bptbl_t *bptbl)
{
    return garray_next_idx(bptbl->ent);
}

int
bptbl_frame_idx(bptbl_t *bptbl)
{
    return bptbl->n_frame;
}

bpidx_t
bptbl_idx(bptbl_t *bptbl, bp_t *bpe)
{
    if (bpe->frame < bptbl->active_fr)
        return bpe - garray_ptr(bptbl->retired, bp_t, 0);
    else
        return bpe - garray_ptr(bptbl->ent, bp_t, garray_base(bptbl->ent))
            + garray_base(bptbl->ent);
}

bp_t *
bptbl_prev(bptbl_t *bptbl, bp_t *ent)
{
    return bptbl_ent(bptbl, ent->bp);
}

int
bptbl_sf(bptbl_t *bptbl, bpidx_t bpidx)
{
    bp_t *ent = bptbl_ent(bptbl, bpidx);
    bp_t *prev;

    if (ent == NULL)
        return -1;
    prev = bptbl_prev(bptbl, ent);
    if (prev == NULL)
        return 0;
    else
        return prev->frame + 1;
}

int
bptbl_ef_count(bptbl_t *bptbl, int frame_idx)
{
    return bptbl_ef_idx(bptbl, frame_idx + 1)
        - bptbl_ef_idx(bptbl, frame_idx);
}

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

bp_t *
bptbl_enter(bptbl_t *bptbl, int32 w, int32 path, int32 score, int rc)
{
    int32 i, rcsize, *bss;
    bp_t be, *bpe;

    /* This might happen if recognition fails. */
    if (bptbl_end_idx(bptbl) == NO_BP) {
        E_ERROR("No entries in backpointer table!");
        return NULL;
    }

    /* Expand the bss table if necessary. */
    if (bptbl->bss_head >= bptbl->bscore_stack_size
        - bin_mdef_n_ciphone(bptbl->d2p->mdef)) {
        bptbl->bscore_stack_size *= 2;
        bptbl->bscore_stack = ckd_realloc(bptbl->bscore_stack,
                                          bptbl->bscore_stack_size
                                          * sizeof(*bptbl->bscore_stack));
        E_INFO("Resized score stack to %d entries (%d KiB)\n",
               bptbl->bscore_stack_size,
               bptbl->bscore_stack_size * sizeof(*bptbl->bscore_stack) / 1024);
    }

    /* Append a new backpointer. */
    memset(&be, 0, sizeof(be));
    be.wid = w;
    be.frame = bptbl->n_frame - 1;
    be.bp = path;
    be.score = score;
    be.s_idx = bptbl->bss_head;
    be.valid = TRUE;
    be.last_phone = dict_last_phone(bptbl->d2p->dict, w);
    bpe = garray_append(bptbl->ent, &be);
    /* Set up its LM state (NOTE: should use the pointer returned by garray_append?). */
    bptbl_fake_lmstate(bptbl, bptbl_end_idx(bptbl) - 1);

    /* DICT2PID */
    /* Get diphone ID for final phone and number of ssids corresponding to it. */
    rcsize = bptbl_rcsize(bptbl, bpe);
    /* Allocate some space on the bptbl->bscore_stack for all of these triphones. */
    for (i = rcsize, bss = bptbl->bscore_stack + bptbl->bss_head; i > 0; --i, bss++)
        *bss = WORST_SCORE;
    bptbl->bscore_stack[bptbl->bss_head + rc] = score;
    bptbl->bss_head += rcsize;

    E_DEBUG(3,("Entered bp %d sf %d ef %d active_fr %d\n", bptbl_end_idx(bptbl) - 1,
               bptbl_sf(bptbl, bptbl_end_idx(bptbl) - 1),
               bptbl->n_frame - 1, bptbl->active_fr));
    assert(bptbl_sf(bptbl, bptbl_end_idx(bptbl) - 1) >= bptbl->active_fr);
    return bpe;
}

void
bptbl_fake_lmstate(bptbl_t *bptbl, int32 bp)
{
    bp_t *ent, *prev;

    assert(bp != NO_BP);
    ent = bptbl_ent(bptbl, bp);
    prev = bptbl_prev(bptbl, ent);
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

char *
bptbl_hyp(bptbl_t *bptbl, int32 *out_score, int32 finish_wid)
{
    bp_t *exit, *bpe;
    char *c, *hyp_str;
    size_t len;

    /* Look for </s> in the last frame. */
    if ((exit = bptbl_find_exit(bptbl, finish_wid)) == NULL) {
        /* If not found then take the best scoring word exit (with warning). */
        exit = bptbl_find_exit(bptbl, BAD_S3WID);
        if (exit == NULL) {
            E_ERROR("No word exits in last frame: recognition failure?\n");
            return NULL;
        }
        E_WARN("No %s found in last frame, using %s instead\n",
               dict_wordstr(bptbl->d2p->dict, finish_wid),
               dict_wordstr(bptbl->d2p->dict, exit->wid));
    }
    if (out_score)
        *out_score = exit->score;

    bpe = exit;
    len = 0;
    while (bpe != NULL) {
        assert(bpe->valid);
        if (dict_real_word(bptbl->d2p->dict, bpe->wid))
            len += strlen(dict_basestr(bptbl->d2p->dict, bpe->wid)) + 1;
        bpe = bptbl_prev(bptbl, bpe);
    }

    if (len == 0)
	return NULL;
    hyp_str = ckd_calloc(1, len);

    bpe = exit;
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
        bpe = bptbl_prev(bptbl, bpe);
    }

    return hyp_str;
}

static void
bptbl_bp2itor(ps_seg_t *seg, int bp)
{
    bptbl_seg_t *bseg = (bptbl_seg_t *)seg;
    bp_t *be, *pbe;

    be = bptbl_ent(bseg->bptbl, bp);
    pbe = bptbl_ent(bseg->bptbl, be->bp);
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
    bp_t *exit, *bpe;
    int cur;

    /* Look for </s> in the last frame. */
    if ((exit = bptbl_find_exit(bptbl, finish_wid)) == NULL) {
        /* If not found then take the best scoring word exit (with warning). */
        exit = bptbl_find_exit(bptbl, BAD_S3WID);
        if (exit == NULL) {
            E_ERROR("No word exits in last frame: recognition failure?\n");
            return NULL;
        }
        E_WARN("No %s found in last frame, using %d instead\n",
               dict_wordstr(bptbl->d2p->dict, finish_wid),
               dict_wordstr(bptbl->d2p->dict, exit->wid));
    }
    if (out_score)
        *out_score = exit->score;

    /* Calling this an "iterator" is a bit of a misnomer since we have
     * to get the entire backtrace in order to produce it.  On the
     * other hand, all we actually need is the bptbl IDs, and we can
     * allocate a fixed-size array of them. */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->base.vt = &bptbl_segfuncs;
    itor->base.lwf = 1.0;
    itor->bptbl = bptbl;
    itor->n_bpidx = 0;
    bpe = exit;
    while (bpe != NULL) {
        ++itor->n_bpidx;
        bpe = bptbl_prev(bptbl, bpe);
    }
    if (itor->n_bpidx == 0) {
        ckd_free(itor);
        return NULL;
    }
    itor->bpidx = ckd_calloc(itor->n_bpidx, sizeof(*itor->bpidx));
    cur = itor->n_bpidx - 1;
    bpe = exit;
    while (bpe != NULL) {
        itor->bpidx[cur] = bptbl_idx(bptbl, bpe);
        bpe = bptbl_prev(bptbl, bpe);
        --cur;
    }

    /* Fill in relevant fields for first element. */
    bptbl_bp2itor((ps_seg_t *)itor, itor->bpidx[0]);

    return (ps_seg_t *)itor;
}

