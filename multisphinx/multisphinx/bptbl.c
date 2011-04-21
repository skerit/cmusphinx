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

#include <multisphinx/bptbl.h>
#include <multisphinx/search.h>
#include "search_internal.h"
#include "hmm.h"

static int bptbl_rcsize(bptbl_t *bptbl, bp_t *be);

#if 0
#undef E_DEBUG
#define E_DEBUG(level,x) E_INFO x
#undef E_DEBUGCONT
#define E_DEBUGCONT(level,x) E_INFOCONT x
#endif

static void
bptbl_lock(bptbl_t *bptbl)
{
    ptmr_start(&bptbl->t_bptbl);
}

static void
bptbl_unlock(bptbl_t *bptbl)
{
    ptmr_stop(&bptbl->t_bptbl);
}

bptbl_t *
bptbl_init(char const *name, dict2pid_t *d2p, int n_alloc, int n_frame_alloc)
{
    bptbl_t *bptbl = ckd_calloc(1, sizeof(*bptbl));

    bptbl->refcount = 1;
    bptbl->name = ckd_salloc(name);
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
    bptbl->rc = garray_init(0, sizeof(rcdelta_t));
    garray_reserve(bptbl->rc, n_alloc * 20); /* 20 = guess at average number of rcs/word */
    ptmr_init(&bptbl->t_bptbl);

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

    E_INFO("%s: TOTAL bptbl %f CPU %f wall\n", bptbl->name,
           bptbl->t_bptbl.t_tot_cpu, bptbl->t_bptbl.t_tot_elapsed);
    dict2pid_free(bptbl->d2p);
    garray_free(bptbl->ent);
    garray_free(bptbl->retired);
    garray_free(bptbl->permute);
    garray_free(bptbl->ef_idx);
    garray_free(bptbl->rc);
    bitvec_free(bptbl->valid_fr);
    ckd_free(bptbl->name);
    ckd_free(bptbl);
    return 0;
}

void
bptbl_reset(bptbl_t *bptbl)
{
    bitvec_clear_all(bptbl->valid_fr, bptbl->n_frame_alloc);
    garray_reset(bptbl->ent);
    garray_reset(bptbl->permute);
    garray_reset(bptbl->ef_idx);
    garray_reset(bptbl->retired);
    garray_reset(bptbl->rc);
    bptbl->dest_s_idx = 0;
    bptbl->n_frame = 0;
    bptbl->oldest_bp = NO_BP;
    ptmr_reset(&bptbl->t_bptbl);
}

void
bptbl_dump(bptbl_t *bptbl)
{
    int i;

    E_INFO("%s: retired backpointers (%d entries, oldest active %d):\n",
           bptbl->name, bptbl_retired_idx(bptbl), bptbl->oldest_bp);
    for (i = garray_base(bptbl->retired);
         i < bptbl_retired_idx(bptbl); ++i) {
        bp_t *ent = garray_ptr(bptbl->retired, bp_t, i);
        assert(ent->valid);
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d real_wid %-5d prev_real_wid %-5d\n",
                    i, dict_wordstr(bptbl->d2p->dict, ent->wid),
                    bptbl_sf(bptbl, i),
                    ent->frame,
                    ent->score,
                    ent->bp,
                    ent->real_wid,
                    ent->prev_real_wid);
    }
    E_INFO("%s: active backpointers (%d entries starting at %d):\n",
           bptbl->name, bptbl_end_idx(bptbl) - bptbl_active_idx(bptbl),
           bptbl_active_idx(bptbl));
    for (i = bptbl_active_idx(bptbl); i < bptbl_end_idx(bptbl); ++i) {
        bp_t *ent = garray_ptr(bptbl->ent, bp_t, i);
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

#ifdef __GNUC__
#define ATTRIBUTE_UNUSED __attribute__((unused))
#else
#define ATTRIBUTE_UNUSED
#endif

/**
 * Mark all active entries in the backpointer table.
 */
static int ATTRIBUTE_UNUSED
bptbl_mark_all(bptbl_t *bptbl, int ef, int cf)
{
    int i, j;

    assert(ef > bptbl_active_frame(bptbl));

    for (i = bptbl_ef_idx(bptbl, ef);
         i < bptbl_ef_idx(bptbl, cf); ++i) {
        bp_t *ent;
        ent = bptbl_ent(bptbl, i);
        ent->valid = TRUE;
    }

    E_DEBUG(2,("Marked bps %d to %d\n",
               bptbl_ef_idx(bptbl, bptbl_active_frame(bptbl)),
               bptbl_ef_idx(bptbl, ef)));
    for (j = 0, i = bptbl_ef_idx(bptbl, bptbl_active_frame(bptbl));
         i < bptbl_ef_idx(bptbl, ef); ++i) {
        bp_t *ent = garray_ptr(bptbl->ent, bp_t, i);
        assert(ent != NULL);
        if (ent->valid)
            ++j;
        E_DEBUG(3,("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                   i, dict_wordstr(bptbl->d2p->dict, ent->wid),
                   bptbl_sf(bptbl, i),
                   ent->frame,
                   ent->score,
                   ent->bp));
    }
    return j;
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

    assert(ef > bptbl_active_frame(bptbl));

    /* Invalidate all active backpointer entries up to ef. */
    E_DEBUG(2,("Invalidating backpointers from %d to %d (%d to %d)\n",
               bptbl_active_frame(bptbl), ef,
               bptbl_ef_idx(bptbl, bptbl_active_frame(bptbl)),
               bptbl_ef_idx(bptbl, ef)));
    for (i = bptbl_ef_idx(bptbl, bptbl_active_frame(bptbl));
         i < bptbl_ef_idx(bptbl, ef); ++i) {
        E_DEBUG(5,("Invalidate bp %d\n", i));
        garray_ent(bptbl->ent, bp_t, i).valid = FALSE;
    }

    /* Now re-activate all ones backwards reachable from the search graph. */
    E_DEBUG(2,("Finding coaccessible frames from backpointers from %d to %d (%d to %d)\n",
               ef, cf,
               bptbl_ef_idx(bptbl, ef), bptbl_ef_idx(bptbl, cf)));
    /* Mark everything immediately reachable from (ef..cf) */
    bitvec_clear_all(bptbl->valid_fr, cf - bptbl_active_frame(bptbl));
    n_active_fr = 0;
    /* NOTE: This for statement can be sped up at the cost of being
     * less obvious. */
    for (i = bptbl_ef_idx(bptbl, ef);
         i < bptbl_ef_idx(bptbl, cf); ++i) {
        bp_t *ent, *prev;
        ent = bptbl_ent(bptbl, i);
        /* This one could be retired. */
        prev = bptbl_ent(bptbl, ent->bp);
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
                for (j = bptbl_ef_idx(bptbl, i);
                     j < bptbl_ef_idx(bptbl, i + 1); ++j) {
                    bp_t *ent = bptbl_ent(bptbl, j);
                    bp_t *prev = bptbl_ent(bptbl, ent->bp);
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
    for (j = 0, i = bptbl_ef_idx(bptbl, bptbl_active_frame(bptbl));
         i < bptbl_ef_idx(bptbl, ef); ++i) {
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
            assert(bptbl_sf(bptbl, i) <= bpe->frame);
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
    if (bptbl_ef_idx(bptbl, bptbl_active_frame(bptbl))
        == bptbl_ef_idx(bptbl, next_active_fr)) {
        bptbl_update_active(bptbl, next_active_fr);
        return;
    }
    E_DEBUG(2,("GC from frame %d to %d\n", bptbl_active_frame(bptbl),
               next_active_fr));
    /* Expand the permutation table if necessary. */
    garray_expand_to(bptbl->permute, bptbl_ef_idx(bptbl, next_active_fr));
    garray_set_base(bptbl->permute, bptbl_active_idx(bptbl));
    /* Mark, compact, snap pointers. */
    n_retired = bptbl_mark(bptbl, next_active_fr, frame_idx);
    E_DEBUG(2,("About to retire %d bps\n", n_retired));
    first_retired_bp = bptbl_retired_idx(bptbl);
    bptbl_retire(bptbl, n_retired,
                 bptbl_ef_idx(bptbl, next_active_fr));
    bptbl_remap(bptbl, first_retired_bp,
                bptbl_ef_idx(bptbl, next_active_fr),
                bptbl_ef_idx(bptbl, next_active_fr));
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
}

int
bptbl_push_frame(bptbl_t *bptbl, int oldest_bp)
{
    int frame_idx = bptbl->n_frame;

    bptbl_lock(bptbl);
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
    bptbl_unlock(bptbl);
    return frame_idx;
}

int
bptbl_commit(bptbl_t *bptbl)
{
    bpidx_t src, dest, eidx;
    int frame_idx, dest_s_idx;

    /* This is the frame we're working on and its bps. */
    bptbl_lock(bptbl);
    frame_idx = bptbl->n_frame - 1;
    dest = bptbl_ef_idx(bptbl, frame_idx);
    dest_s_idx = garray_ent(bptbl->ent, bp_t, dest).s_idx;
    eidx = bptbl_end_idx(bptbl);
    E_DEBUG(4,("compacting %d bps\n", eidx - dest));
    if (eidx == dest) {
        bptbl_unlock(bptbl);
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
                eidx - bptbl_ef_idx(bptbl, frame_idx)));

    /* Truncate active arrays. */
    garray_pop_from(bptbl->rc, dest_s_idx);
    garray_pop_from(bptbl->ent, dest);
    bptbl_unlock(bptbl);

    return dest - src;
}

int
bptbl_is_final(bptbl_t *bptbl)
{
    return (bptbl_end_idx(bptbl) == bptbl_active_idx(bptbl));
}

int
bptbl_finalize(bptbl_t *bptbl)
{
    int n_retired, first_retired_bp;

    bptbl_lock(bptbl);
    E_DEBUG(2,("Final GC from frame %d to %d\n",
               bptbl_active_frame(bptbl), bptbl->n_frame));
    /* If there is nothing to GC then finish up. */
    if (bptbl_is_final(bptbl)) {
        bptbl_unlock(bptbl);
        return 0;
    }
    /* Expand the permutation table if necessary (probably). */
    garray_expand_to(bptbl->permute, bptbl_end_idx(bptbl));
    garray_set_base(bptbl->permute, bptbl_active_idx(bptbl));
    /* Mark and GC everything from the last frame. */
    n_retired = bptbl_mark(bptbl, bptbl->n_frame - 1, bptbl->n_frame);
    /* Include the last frame in the retired count. */
    n_retired += bptbl_ef_count(bptbl, bptbl->n_frame - 1);
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
    E_INFO("%s: allocated %d active and %d retired entries (%d + %d KiB)\n",
           bptbl->name,
           garray_alloc_size(bptbl->ent),
           garray_alloc_size(bptbl->retired),
           garray_alloc_size(bptbl->ent) * sizeof(bp_t) / 1024,
           garray_alloc_size(bptbl->retired) * sizeof(bp_t) / 1024);
    E_INFO("%s: allocated %d right context deltas (%d KiB)\n",
           bptbl->name,
           garray_alloc_size(bptbl->rc),
           garray_alloc_size(bptbl->rc) * sizeof(rcdelta_t) / 1024);
    E_INFO("%s: allocated %d permutation entries and %d end frame entries\n",
           bptbl->name,
           garray_alloc_size(bptbl->permute), garray_alloc_size(bptbl->ef_idx));
    bptbl_unlock(bptbl);
    E_INFO("%s: bptbl %f CPU %f wall %f xRT\n",
           bptbl->name,
           bptbl->t_bptbl.t_cpu, bptbl->t_bptbl.t_elapsed,
           bptbl->t_bptbl.t_elapsed / bptbl->n_frame * 100);
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

    bptbl_lock(bptbl);
    if (first_idx > bptbl_retired_idx(bptbl)) {
        E_DEBUG(2,("%d outside retired, releasing up to %d\n",
                   first_idx, bptbl_retired_idx(bptbl)));
        first_idx = bptbl_retired_idx(bptbl);
    }

    base_idx = garray_base(bptbl->retired);
    E_DEBUG(2,("Releasing bptrs from %d to %d:\n",
               base_idx, first_idx));
    if (first_idx < base_idx) {
        bptbl_unlock(bptbl);
        return 0;
    }
#if 0 /* For debugging purposes... */
    int i;
    for (i = base_idx; i < first_idx; ++i) {
        bp_t *ent = garray_ptr(bptbl->retired, bp_t, i);
        assert(ent->valid);
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                    i, dict_wordstr(bptbl->d2p->dict, ent->wid),
                    bptbl_sf(bptbl, i),
                    ent->frame,
                    ent->score,
                    ent->bp);
    }
#endif

    ent = garray_ptr(bptbl->retired, bp_t, first_idx);
    garray_shift_from(bptbl->rc, ent->s_idx);
    garray_set_base(bptbl->rc, ent->s_idx);
    garray_shift_from(bptbl->retired, first_idx);
    garray_set_base(bptbl->retired, first_idx);

    bptbl_unlock(bptbl);
    return first_idx - base_idx;
}

bpidx_t
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
        start = end = bptbl_retired_idx(bptbl) - 1;
        ef = bptbl_ent(bptbl, start)->frame;
        while (start >= first_retired) {
            if (bptbl_ent(bptbl, start)->frame != ef)
                break;
            --start;
        }
    }
    else {
        bpidx_t first_ent = garray_base(bptbl->ent);
        /* Not final, so it's in ent. */
        start = end = garray_next_idx(bptbl->ent) - 1;
        ef = bptbl_ent(bptbl, start)->frame;
        while (start >= first_ent) {
            if (bptbl_ent(bptbl, start)->frame != ef)
                break;
            --start;
        }
    }
    ++start;
    if (ef != bptbl->n_frame - 1) {
#if 0 /* shut up */
        E_WARN("No exits in final frame %d, using frame %d instead\n",
               bptbl->n_frame - 1, ef);
#endif
    }
    best_score = WORST_SCORE;
    best = NO_BP;
    while (start <= end) {
        bp_t *ent = bptbl_ent(bptbl, start);
        if (ent->score BETTER_THAN best_score
            && (wid == BAD_S3WID || ent->wid == wid)) {
            best = start;
            best_score = ent->score;
        }
        ++start;
    }
    return best;
}

int32
bptbl_ef_idx(bptbl_t *bptbl, int frame_idx)
{
    if (frame_idx < bptbl_active_frame(bptbl))
        return 0;
    else if (frame_idx >= bptbl->n_frame)
        return bptbl_end_idx(bptbl);
    else
        return garray_ent(bptbl->ef_idx, bpidx_t, frame_idx);
}

bp_t *
bptbl_ent(bptbl_t *bptbl, bpidx_t bpidx)
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

    ent = bptbl_ent(bptbl, bpidx);
    if (ent == NULL) {
        return -1;
    }
    memcpy(out_bp, ent, sizeof(*out_bp));

    return 0;
}

int
bptbl_set_bp(bptbl_t *bptbl, bpidx_t bpidx, bp_t const *bp)
{
    bp_t *ent;

    if (bpidx == NO_BP)
        return -1;

    ent = bptbl_ent(bptbl, bpidx);
    if (ent == NULL) {
        return -1;
    }
    memcpy(ent, bp, sizeof(*ent));

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

    if (bptbl->oldest_bp == NO_BP)
        sf = 0;
    else {
        bp_t *ent;
        ent = bptbl_ent(bptbl, bptbl->oldest_bp);
        sf = ent->frame + 1;
    }
    bptbl_unlock(bptbl);
    return sf;
}

int
bptbl_sf(bptbl_t *bptbl, bpidx_t bpidx)
{
    bp_t *ent;
    bp_t *prev;

    ent = bptbl_ent(bptbl, bpidx);
    if (ent == NULL)
        return -1;
    else {
        prev = bptbl_ent(bptbl, ent->bp);
        if (prev == NULL)
            return 0;
        else
            return prev->frame + 1;
    }
}

int
bptbl_ef_count(bptbl_t *bptbl, int frame_idx)
{
    bpidx_t start, end;
    start = bptbl_ef_idx(bptbl, frame_idx);
    end = bptbl_ef_idx(bptbl, frame_idx + 1);
    return end - start;
}

void
bptbl_set_rcscore(bptbl_t *bptbl, bpidx_t bpidx, int rc, int32 score)
{
    bp_t *bpe;
    int delta;

    bpe = bptbl_ent(bptbl, bpidx);
    /* No rc scores to set! */
    if (dict_is_single_phone(bptbl->d2p->dict, bpe->wid))
        return;
    assert(score <= bpe->score);
    delta = bpe->score - score;
    if (score == WORST_SCORE || delta >= NO_RC)
        garray_ent(bptbl->rc, rcdelta_t, bpe->s_idx + rc) = NO_RC;
    else
        garray_ent(bptbl->rc, rcdelta_t, bpe->s_idx + rc) = delta;
}

int
bptbl_get_rcscores(bptbl_t *bptbl, bpidx_t bpidx, int32 *out_rcscores)
{
    bp_t *bpe;
    int rcsize;

    bpe = bptbl_ent(bptbl, bpidx);
    rcsize = bptbl_rcsize(bptbl, bpe);
    if (rcsize == 0) {
        assert(dict_is_single_phone(bptbl->d2p->dict, bpe->wid));
        out_rcscores[0] = bpe->score;
        return 1;
    }
    else {
        int i;
        assert(bpe->s_idx < garray_next_idx(bptbl->rc));
        for (i = 0; i < rcsize; ++i) {
            int delta = garray_ent(bptbl->rc, rcdelta_t, bpe->s_idx + i);
            if (delta == NO_RC)
                out_rcscores[i] = WORST_SCORE;
            else
                out_rcscores[i] = bpe->score - delta;
        }
        return rcsize;
    }
}

int
bptbl_get_rcdeltas(bptbl_t *bptbl, bpidx_t bpidx, rcdelta_t *out_rcdeltas)
{
    bp_t *bpe;
    int rcsize;

    bpe = bptbl_ent(bptbl, bpidx);
    rcsize = bptbl_rcsize(bptbl, bpe);
    if (rcsize == 0) {
        assert(dict_is_single_phone(bptbl->d2p->dict, bpe->wid));
        out_rcdeltas[0] = 0;
        return 1;
    }
    else {
        memcpy(out_rcdeltas,
               garray_ptr(bptbl->rc, rcdelta_t, bpe->s_idx),
               rcsize * sizeof(rcdelta_t));
        return rcsize;
    }
}

/**
 * Get the number of right context entries for a backpointer.
 *
 * This is only used internally.  It returns the correct rcsize of
 * zero rather than a "cooked" value of one (see bptbl_get_rcscores())
 * in the case where there are no right contexts.
 */
static int
bptbl_rcsize(bptbl_t *bptbl, bp_t *be)
{
    int rcsize;

    if (dict_is_single_phone(bptbl->d2p->dict, be->wid)) {
        be->last2_phone = -1;
        rcsize = 0;
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

    prev = bptbl_ent(bptbl, ent->bp);
    /* Propagate lm state for fillers, rotate it for words. */
    if (dict_filler_word(bptbl->d2p->dict, ent->wid)) {
        if (prev != NULL) {
            ent->real_wid = prev->real_wid;
            ent->prev_real_wid = prev->prev_real_wid;
        }
        else {
            ent->real_wid = dict_basewid(bptbl->d2p->dict, ent->wid);
            ent->prev_real_wid = BAD_S3WID;
        }
    }
    else {
        ent->real_wid = dict_basewid(bptbl->d2p->dict, ent->wid);
        if (prev != NULL)
            ent->prev_real_wid = prev->real_wid;
        else
            ent->prev_real_wid = BAD_S3WID;
    }
}

int32
bptbl_fake_lmscore(bptbl_t *bptbl, ngram_model_t *lm,
                   bpidx_t bp, int *out_n_used)
{
    bp_t *prev, *ent;

    ent = bptbl_ent(bptbl, bp);
    assert(ent != NULL);
    prev = bptbl_ent(bptbl, ent->bp);

    /* Start word has lscr = 0 */
    if (prev == NULL)
        return 0;
    /* Filler/silence words have search-dependent scores. */
    else if (dict_filler_word(bptbl->d2p->dict, ent->wid)
             || ent->wid == dict_startwid(bptbl->d2p->dict))
        return 0;
    else {
        return ngram_tg_score(lm, ent->real_wid,
                              prev->real_wid,
                              prev->prev_real_wid,
                              out_n_used) >> SENSCR_SHIFT;
    }
}

bpidx_t
bptbl_enter(bptbl_t *bptbl, int32 w, int32 path, int32 score, int rc)
{
    int32 i, rcsize;
    size_t bss_head;
    bpidx_t bpidx;
    bp_t be, *bpe;
    rcdelta_t *bss;

    /* This might happen if recognition fails. */
    if (bptbl_end_idx(bptbl) == NO_BP) {
        E_ERROR("No entries in backpointer table!");
        return NO_BP;
    }

    bptbl_lock(bptbl);
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
    if (rcsize) {
        bss_head = be.s_idx;
        garray_expand_to(bptbl->rc, bss_head + rcsize);
        bss = garray_ptr(bptbl->rc, rcdelta_t, bss_head);
        for (i = 0; i < rcsize; ++i)
            *bss++ = NO_RC;
        garray_ent(bptbl->rc, rcdelta_t, bss_head + rc) = 0;
    }
    E_DEBUG(3,("Entered bp %d sf %d ef %d s_idx %d active_fr %d\n",
               bptbl_end_idx(bptbl) - 1,
               bptbl_sf(bptbl, bptbl_end_idx(bptbl) - 1),
               bptbl->n_frame - 1, bss_head, bptbl_active_frame(bptbl)));
    assert(bptbl_sf(bptbl, bptbl_end_idx(bptbl) - 1)
           >= bptbl_active_frame(bptbl));
    bptbl_unlock(bptbl);
    return bpidx;
}

void
bptbl_update_bp(bptbl_t *bptbl, int32 bp, int rc,
                bpidx_t new_prev, int32 new_score)
{
    bp_t *ent;
    int rcsize;

    assert(bp != NO_BP);
    ent = bptbl_ent(bptbl, bp);
    assert(new_score > ent->score);
    rcsize = bptbl_rcsize(bptbl, ent);
    if (rcsize != 0) {
        int i, delta;
        assert(ent->s_idx < garray_next_idx(bptbl->rc));
        delta = new_score - ent->score;
        if (delta > 0) {
            for (i = 0; i < rcsize; ++i) {
                int cur_score = garray_ent(bptbl->rc, rcdelta_t, ent->s_idx + i);
                if (cur_score == NO_RC)
                    continue;
                else if (cur_score + delta >= NO_RC) {
                    E_WARN("rc score overflow in bp %d rc %d: %d + %d\n",
                           bp, i, cur_score, delta);
                    garray_ent(bptbl->rc, rcdelta_t, ent->s_idx + i) = NO_RC;
                }
                else
                    garray_ent(bptbl->rc, rcdelta_t, ent->s_idx + i) += delta;
            }
        }
    }
    ent->bp = new_prev;
    ent->score = new_score;
    bptbl_fake_lmstate_internal(bptbl, ent);
}

char *
bptbl_backtrace(bptbl_t *bptbl, bpidx_t bp)
{
    char *c, *hyp_str;
    size_t len;
    bp_t *bpe;

    bpe = bptbl_ent(bptbl, bp);
    len = 0;

    while (bpe != NULL) {
        assert(bpe->valid);
        if (dict_real_word(bptbl->d2p->dict, bpe->wid))
            len += strlen(dict_basestr(bptbl->d2p->dict, bpe->wid)) + 1;
        bpe = bptbl_ent(bptbl, bpe->bp);
    }

    if (len == 0) {
        return NULL;
    }
    hyp_str = ckd_calloc(1, len);

    bpe = bptbl_ent(bptbl, bp);
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
        bpe = bptbl_ent(bptbl, bpe->bp);
    }

    return hyp_str;
}

char *
bptbl_hyp(bptbl_t *bptbl, int32 *out_score, int32 finish_wid)
{
    bpidx_t exit;
    bp_t *bpe;

    /* Look for </s> in the last frame. */
    if ((exit = bptbl_find_exit(bptbl, finish_wid)) == NO_BP) {
        /* If not found then take the best scoring word exit (with warning). */
        exit = bptbl_find_exit(bptbl, BAD_S3WID);
        if (exit == NO_BP) {
            E_ERROR("No word exits in last frame: recognition failure?\n");
            return NULL;
        }
#if 0 /* shut up */
        E_WARN("No %s found in last frame, using %s instead\n",
               dict_wordstr(bptbl->d2p->dict, finish_wid),
               dict_wordstr(bptbl->d2p->dict,
                            bptbl_ent(bptbl, exit)->wid));
#endif
    }

    bpe = bptbl_ent(bptbl, exit);
    if (out_score)
        *out_score = bpe->score;
    return bptbl_backtrace(bptbl, exit);
}


/**
 * Segmentation "iterator" for backpointer table results.
 */
typedef struct bptbl_seg_s {
    struct seg_iter_s base;  /**< Base structure. */
    bptbl_t *bptbl; /**< Backpointer table. */
    int32 *bpidx;   /**< Sequence of backpointer IDs. */
    int16 n_bpidx;  /**< Number of backpointer IDs. */
    int16 cur;      /**< Current position in bpidx. */
} bptbl_seg_t;

static void
bptbl_bp2itor(seg_iter_t *seg, int bp)
{
    bptbl_seg_t *bseg = (bptbl_seg_t *)seg;
    bp_t *be, *pbe;

    be = bptbl_ent(bseg->bptbl, bp);
    pbe = bptbl_ent(bseg->bptbl, be->bp);
    seg->wid = be->wid;
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
bptbl_seg_free(seg_iter_t *seg)
{
    bptbl_seg_t *itor = (bptbl_seg_t *)seg;
    
    ckd_free(itor->bpidx);
    ckd_free(itor);
}

static seg_iter_t *
bptbl_seg_next(seg_iter_t *seg)
{
    bptbl_seg_t *itor = (bptbl_seg_t *)seg;

    if (++itor->cur == itor->n_bpidx) {
        bptbl_seg_free(seg);
        return NULL;
    }

    bptbl_bp2itor(seg, itor->bpidx[itor->cur]);
    return seg;
}

static segfuncs_t bptbl_segfuncs = {
    /* seg_next */ bptbl_seg_next,
    /* seg_free */ bptbl_seg_free
};

seg_iter_t *
bptbl_seg_iter(bptbl_t *bptbl, int32 *out_score, int32 finish_wid)
{
    bpidx_t bp;
    bp_t *bpe;

    /* Look for </s> in the last frame. */
    if ((bp = bptbl_find_exit(bptbl, finish_wid)) == NO_BP) {
        /* If not found then take the best scoring word exit (with warning). */
        bp = bptbl_find_exit(bptbl, BAD_S3WID);
        if (bp == NO_BP) {
            E_ERROR("No word exits in last frame: recognition failure?\n");
            return NULL;
        }
        E_WARN("No %s found in last frame, using %s instead\n",
               dict_wordstr(bptbl->d2p->dict, finish_wid),
               dict_wordstr(bptbl->d2p->dict,
                            bptbl_ent(bptbl, bp)->wid));
    }

    bpe = bptbl_ent(bptbl, bp);
    if (out_score)
        *out_score = bpe->score;

    return bptbl_seg_backtrace(bptbl, bp);
}

seg_iter_t *
bptbl_seg_backtrace(bptbl_t *bptbl, bpidx_t bp)
{
    bptbl_seg_t *itor;
    bp_t *bpe;
    int cur;

    bpe = bptbl_ent(bptbl, bp);

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
        bpe = bptbl_ent(bptbl, bpe->bp);
    }
    if (itor->n_bpidx == 0) {
        ckd_free(itor);
        return NULL;
    }
    itor->bpidx = ckd_calloc(itor->n_bpidx, sizeof(*itor->bpidx));
    cur = itor->n_bpidx - 1;
    bpe = bptbl_ent(bptbl, bp);
    while (bpe != NULL) {
        itor->bpidx[cur] = bp;
        bp = bpe->bp;
        bpe = bptbl_ent(bptbl, bp);
        --cur;
    }


    /* Fill in relevant fields for first element. */
    bptbl_bp2itor((seg_iter_t *)itor, itor->bpidx[0]);

    return (seg_iter_t *)itor;
}
