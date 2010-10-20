/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2010 Carnegie Mellon University.  All rights
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

/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#if defined(__ADSPBLACKFIN__)
#elif !defined(_WIN32_WCE)
#include <sys/types.h>
#endif

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

/* SphinxBase headers */
#include <sphinx_config.h>
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/fixpoint.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/bio.h>
#include <sphinxbase/err.h>
#include <sphinxbase/prim_type.h>

/* Local headers */
#include "tied_mgau_common.h"
#include "ptm_mgau.h"
#include "posixwin32.h"

static ps_mgaufuncs_t ptm_mgau_funcs = {
    "ptm",
    &ptm_mgau_frame_eval,      /* frame_eval */
    &ptm_mgau_copy,            /* copy */
    &ptm_mgau_free             /* free */
};

#define COMPUTE_GMM_MAP(_idx)                           \
    diff[_idx] = obs[_idx] - mean[_idx];                \
    sqdiff[_idx] = MFCCMUL(diff[_idx], diff[_idx]);     \
    compl[_idx] = MFCCMUL(sqdiff[_idx], var[_idx]);
#define COMPUTE_GMM_REDUCE(_idx)                \
    d = GMMSUB(d, compl[_idx]);

static void
insertion_sort_topn(ptm_topn_t *topn, int i, int32 d)
{
    ptm_topn_t vtmp;
    int j;

    topn[i].score = d;
    if (i == 0)
        return;
    vtmp = topn[i];
    for (j = i - 1; j >= 0 && d > topn[j].score; j--) {
        topn[j + 1] = topn[j];
    }
    topn[j + 1] = vtmp;
}

static int
eval_topn(ptm_mgau_t *s, int cb, int feat, mfcc_t *z)
{
    ptm_topn_t *topn;
    int i, ceplen;

    topn = s->f->topn[cb][feat];
    ceplen = s->g->featlen[feat];

    for (i = 0; i < s->max_topn; i++) {
        mfcc_t *mean, diff[4], sqdiff[4], compl[4]; /* diff, diff^2, component likelihood */
        mfcc_t *var, d;
        mfcc_t *obs;
        int32 cw, j;

        cw = topn[i].cw;
        mean = s->g->mean[cb][feat][0] + cw * ceplen;
        var = s->g->var[cb][feat][0] + cw * ceplen;
        d = s->g->det[cb][feat][cw];
        obs = z;
        for (j = 0; j < ceplen % 4; ++j) {
            diff[0] = *obs++ - *mean++;
            sqdiff[0] = MFCCMUL(diff[0], diff[0]);
            compl[0] = MFCCMUL(sqdiff[0], *var);
            d = GMMSUB(d, compl[0]);
            ++var;
        }
        /* We could vectorize this but it's unlikely to make much
         * difference as the outer loop here isn't very big. */
        for (;j < ceplen; j += 4) {
            COMPUTE_GMM_MAP(0);
            COMPUTE_GMM_MAP(1);
            COMPUTE_GMM_MAP(2);
            COMPUTE_GMM_MAP(3);
            COMPUTE_GMM_REDUCE(0);
            COMPUTE_GMM_REDUCE(1);
            COMPUTE_GMM_REDUCE(2);
            COMPUTE_GMM_REDUCE(3);
            var += 4;
            obs += 4;
            mean += 4;
        }
        insertion_sort_topn(topn, i, (int32)d);
    }

    return topn[0].score;
}

/* This looks bad, but it actually isn't.  Less than 1% of eval_cb's
 * time is spent doing this. */
static void
insertion_sort_cb(ptm_topn_t **cur, ptm_topn_t *worst, ptm_topn_t *best,
                  int cw, int32 intd)
{
    for (*cur = worst - 1; *cur >= best && intd >= (*cur)->score; --*cur)
        memcpy(*cur + 1, *cur, sizeof(**cur));
    ++*cur;
    (*cur)->cw = cw;
    (*cur)->score = intd;
}

static int
eval_cb(ptm_mgau_t *s, int cb, int feat, mfcc_t *z)
{
    ptm_topn_t *worst, *best, *topn;
    mfcc_t *mean;
    mfcc_t *var, *det, *detP, *detE;
    int32 i, ceplen;

    best = topn = s->f->topn[cb][feat];
    worst = topn + (s->max_topn - 1);
    mean = s->g->mean[cb][feat][0];
    var = s->g->var[cb][feat][0];
    det = s->g->det[cb][feat];
    detE = det + s->g->n_density;
    ceplen = s->g->featlen[feat];

    for (detP = det; detP < detE; ++detP) {
        mfcc_t diff[4], sqdiff[4], compl[4]; /* diff, diff^2, component likelihood */
        mfcc_t d, thresh;
        mfcc_t *obs;
        ptm_topn_t *cur;
        int32 cw, j;

        d = *detP;
        thresh = (mfcc_t) worst->score; /* Avoid int-to-float conversions */
        obs = z;
        cw = detP - det;

        /* Unroll the loop starting with the first dimension(s).  In
         * theory this might be a bit faster if this Gaussian gets
         * "knocked out" by C0. In practice not. */
        for (j = 0; (j < ceplen % 4) && (d >= thresh); ++j) {
            diff[0] = *obs++ - *mean++;
            sqdiff[0] = MFCCMUL(diff[0], diff[0]);
            compl[0] = MFCCMUL(sqdiff[0], *var++);
            d = GMMSUB(d, compl[0]);
        }
        /* Now do 4 dimensions at a time.  You'd think that GCC would
         * vectorize this?  Apparently not.  And it's right, because
         * that won't make this any faster, at least on x86-64. */
        for (; j < ceplen && d >= thresh; j += 4) {
            COMPUTE_GMM_MAP(0);
            COMPUTE_GMM_MAP(1);
            COMPUTE_GMM_MAP(2);
            COMPUTE_GMM_MAP(3);
            COMPUTE_GMM_REDUCE(0);
            COMPUTE_GMM_REDUCE(1);
            COMPUTE_GMM_REDUCE(2);
            COMPUTE_GMM_REDUCE(3);
            var += 4;
            obs += 4;
            mean += 4;
        }
        if (j < ceplen) {
            /* terminated early, so not in topn */
            mean += (ceplen - j);
            var += (ceplen - j);
            continue;
        }
        if (d < thresh)
            continue;
        for (i = 0; i < s->max_topn; i++) {
            /* already there, so don't need to insert */
            if (topn[i].cw == cw)
                break;
        }
        if (i < s->max_topn)
            continue;       /* already there.  Don't insert */
        insertion_sort_cb(&cur, worst, best, cw, (int32)d);
    }

    return best->score;
}

/**
 * Compute top-N densities for active codebooks (and prune)
 */
static int
ptm_mgau_codebook_eval(ptm_mgau_t *s, mfcc_t **z, int frame)
{
    int i, j;

    /* First evaluate top-N from previous frame. */
    for (i = 0; i < s->g->n_mgau; ++i)
        for (j = 0; j < s->g->n_feat; ++j)
            eval_topn(s, i, j, z[j]);

    /* If frame downsampling is in effect, possibly do nothing else. */
    if (frame % s->ds_ratio)
        return 0;

    /* Evaluate remaining codebooks. */
    for (i = 0; i < s->g->n_mgau; ++i) {
        if (bitvec_is_clear(s->f->mgau_active, i))
            continue;
        for (j = 0; j < s->g->n_feat; ++j) {
            eval_cb(s, i, j, z[j]);
        }
    }

    /* Normalize densities to produce "posterior probabilities",
     * i.e. things with a reasonable dynamic range, then scale and
     * clamp them to the acceptable range.  This is actually done
     * solely to ensure that we can use fast_logmath_add().  Note that
     * unless we share the same normalizer across all codebooks for
     * each feature stream we get defective scores (that's why these
     * loops are inside out - doing it per-feature should give us
     * greater precision). */
    for (j = 0; j < s->g->n_feat; ++j) {
        int32 norm = 0x7fffffff;
        for (i = 0; i < s->g->n_mgau; ++i) {
            if (bitvec_is_clear(s->f->mgau_active, i))
                continue;
            if (norm > s->f->topn[i][j][0].score >> SENSCR_SHIFT)
                norm = s->f->topn[i][j][0].score >> SENSCR_SHIFT;
        }
        assert(norm != 0x7fffffff);
        for (i = 0; i < s->g->n_mgau; ++i) {
            int32 k;
            if (bitvec_is_clear(s->f->mgau_active, i))
                continue;
            for (k = 0; k < s->max_topn; ++k) {
                s->f->topn[i][j][k].score >>= SENSCR_SHIFT;
                s->f->topn[i][j][k].score -= norm;
                s->f->topn[i][j][k].score = -s->f->topn[i][j][k].score;
                if (s->f->topn[i][j][k].score > MAX_NEG_ASCR) 
                    s->f->topn[i][j][k].score = MAX_NEG_ASCR;
            }
        }
    }

    return 0;
}

static int
ptm_mgau_calc_cb_active(ptm_mgau_t *s, uint8 *senone_active,
                        int32 n_senone_active, int compallsen)
{
    int i, lastsen;

    if (compallsen) {
        bitvec_set_all(s->f->mgau_active, s->g->n_mgau);
        return 0;
    }
    bitvec_clear_all(s->f->mgau_active, s->g->n_mgau);
    for (lastsen = i = 0; i < n_senone_active; ++i) {
        int sen = senone_active[i] + lastsen;
        int cb = s->s->sen2cb[sen];
        bitvec_set(s->f->mgau_active, cb);
        lastsen = sen;
    }
    E_DEBUG(1, ("Active codebooks:"));
    for (i = 0; i < s->g->n_mgau; ++i) {
        if (bitvec_is_clear(s->f->mgau_active, i))
            continue;
        E_DEBUGCONT(1, (" %d", i));
    }
    E_DEBUGCONT(1, ("\n"));
    return 0;
}

/**
 * Compute senone scores from top-N densities for active codebooks.
 */
static int
ptm_mgau_senone_eval(ptm_mgau_t *s, int16 *senone_scores,
                     uint8 *senone_active, int32 n_senone_active,
                     int compall)
{
    int i, lastsen, bestscore;

    memset(senone_scores, 0, s->n_sen * sizeof(*senone_scores));
    /* FIXME: This is the non-cache-efficient way to do this.  We want
     * to evaluate one codeword at a time but this requires us to have
     * a reverse codebook to senone mapping, which we don't have
     * (yet), since different codebooks have different top-N
     * codewords. */
    if (compall)
        n_senone_active = s->n_sen;
    bestscore = 0x7fffffff;
    for (lastsen = i = 0; i < n_senone_active; ++i) {
        int sen, f, cb;
        int ascore;

        if (compall)
            sen = i;
        else
            sen = senone_active[i] + lastsen;
        lastsen = sen;
        cb = s->s->sen2cb[sen];

        if (bitvec_is_clear(s->f->mgau_active, cb)) {
            int j;
            /* Because senone_active is deltas we can't really "knock
             * out" senones from pruned codebooks, and in any case,
             * it wouldn't make any difference to the search code,
             * which doesn't expect senone_active to change. */
            for (f = 0; f < s->g->n_feat; ++f) {
                for (j = 0; j < s->max_topn; ++j) {
                    s->f->topn[cb][f][j].score = MAX_NEG_ASCR;
                }
            }
        }
        /* For each feature, log-sum codeword scores + mixw to get
         * feature density, then sum (multiply) to get ascore */
        ascore = 0;
        for (f = 0; f < s->g->n_feat; ++f) {
            ptm_topn_t *topn;
            int j, fden = 0;
            topn = s->f->topn[cb][f];
            for (j = 0; j < s->max_topn; ++j) {
                int mixw;
                /* Find mixture weight for this codeword. */
                if (s->s->mixw_cb) {
                    int dcw = s->s->mixw[f][topn[j].cw][sen/2];
                    dcw = (dcw & 1) ? dcw >> 4 : dcw & 0x0f;
                    mixw = s->s->mixw_cb[dcw];
                }
                else {
                    mixw = s->s->mixw[f][topn[j].cw][sen];
                }
                if (j == 0)
                    fden = mixw + topn[j].score;
                else
                    fden = fast_logmath_add(s->lmath_8b, fden,
                                       mixw + topn[j].score);
                E_DEBUG(3, ("fden[%d][%d] l+= %d + %d = %d\n",
                            sen, f, mixw, topn[j].score, fden));
            }
            ascore += fden;
        }
        if (ascore < bestscore) bestscore = ascore;
        senone_scores[sen] = ascore;
    }
    /* Normalize the scores again (finishing the job we started above
     * in ptm_mgau_codebook_eval...) */
    for (i = 0; i < s->n_sen; ++i) {
        senone_scores[i] -= bestscore;
    }

    return 0;
}

/**
 * Compute senone scores for the active senones.
 */
int32
ptm_mgau_frame_eval(ps_mgau_t *ps,
                    int16 *senone_scores,
                    uint8 *senone_active,
                    int32 n_senone_active,
                    mfcc_t ** featbuf, int32 frame,
                    int32 compallsen)
{
    ptm_mgau_t *s = (ptm_mgau_t *)ps;
    int fast_eval_idx;

    /* Find the appropriate frame in the rotating history buffer
     * corresponding to the requested input frame.  No bounds checking
     * is done here, which just means you'll get semi-random crap if
     * you request a frame in the future or one that's too far in the
     * past.  Since the history buffer is just used for fast match
     * that might not be fatal. */
    fast_eval_idx = frame % s->n_fast_hist;
    s->f = s->hist + fast_eval_idx;
    /* Compute the top-N codewords for every codebook, unless this
     * is a past frame, in which case we already have them (we
     * hope!) */
    if (frame >= ps_mgau_base(ps)->frame_idx) {
        ptm_fast_eval_t *lastf;
        /* Get the previous frame's top-N information (on the
         * first frame of the input this is just all WORST_DIST,
         * no harm in that) */
        if (fast_eval_idx == 0)
            lastf = s->hist + s->n_fast_hist - 1;
        else
            lastf = s->hist + fast_eval_idx - 1;
        /* Copy in initial top-N info */
        memcpy(s->f->topn[0][0], lastf->topn[0][0],
               s->g->n_mgau * s->g->n_feat * s->max_topn * sizeof(ptm_topn_t));
        /* Generate initial active codebook list (this might not be
         * necessary) */
        ptm_mgau_calc_cb_active(s, senone_active, n_senone_active, compallsen);
        /* Now evaluate top-N, prune, and evaluate remaining codebooks. */
        ptm_mgau_codebook_eval(s, featbuf, frame);
    }
    /* Evaluate intersection of active senones and active codebooks. */
    ptm_mgau_senone_eval(s, senone_scores, senone_active,
                         n_senone_active, compallsen);

    return 0;
}

ps_mgau_t *
ptm_mgau_init(acmod_t *acmod)
{
    ptm_mgau_t *s;
    ps_mgau_t *ps;
    char const *sendump_path;
    int i;

    s = ckd_calloc(1, sizeof(*s));
    s->config = acmod->config;

    s->lmath = logmath_retain(acmod->lmath);
    /* Log-add table. */
    s->lmath_8b = logmath_init(logmath_get_base(acmod->lmath), SENSCR_SHIFT, TRUE);
    if (s->lmath_8b == NULL)
        goto error_out;
    /* Ensure that it is only 8 bits wide so that fast_logmath_add() works. */
    if (logmath_get_width(s->lmath_8b) != 1) {
        E_ERROR("Log base %f is too small to represent add table in 8 bits\n",
                logmath_get_base(s->lmath_8b));
        goto error_out;
    }

    /* Read means and variances. */
    if ((s->g = gauden_init(cmd_ln_str_r(s->config, "-mean"),
                            cmd_ln_str_r(s->config, "-var"),
                            cmd_ln_float32_r(s->config, "-varfloor"),
                            s->lmath)) == NULL)
        goto error_out;
    /* We only support 256 codebooks or less (like 640k or 2GB, this
     * should be enough for anyone) */
    if (s->g->n_mgau > 256) {
        E_ERROR("Number of codebooks exceeds 256: %d\n", s->g->n_mgau);
        goto error_out;
    }
    /* Verify n_feat and veclen, against acmod. */
    if (s->g->n_feat != feat_dimension1(acmod->fcb)) {
        E_ERROR("Number of streams does not match: %d != %d\n",
                s->g->n_feat, feat_dimension(acmod->fcb));
        goto error_out;
    }
    for (i = 0; i < s->g->n_feat; ++i) {
        if (s->g->featlen[i] != feat_dimension2(acmod->fcb, i)) {
            E_ERROR("Dimension of stream %d does not match: %d != %d\n",
                    s->g->featlen[i], feat_dimension2(acmod->fcb, i));
            goto error_out;
        }
    }
    /* Read mixture weights. */
    if ((sendump_path = cmd_ln_str_r(s->config, "-sendump"))) {
        s->s = sendump_read_sendump(s->config, s->lmath_8b,
                                    s->g, acmod->mdef, sendump_path);
        if (s->s == NULL)
            goto error_out;
    }
    else {
        s->s = sendump_read_mixw(s->config, s->lmath_8b,
                                 s->g, acmod->mdef, sendump_path);
        if (s->s == NULL)
            goto error_out;
    }
    s->n_sen = bin_mdef_n_sen(acmod->mdef);
    s->ds_ratio = cmd_ln_int32_r(s->config, "-ds");
    s->max_topn = cmd_ln_int32_r(s->config, "-topn");
    E_INFO("Maximum top-N: %d\n", s->max_topn);

    /* Assume mapping of senones to their base phones, though this
     * will become more flexible in the future. */
    for (i = 0; i < s->n_sen; ++i)
        s->s->sen2cb[i] = bin_mdef_sen2cimap(acmod->mdef, i);

    /* Allocate fast-match history buffers.  We need enough for the
     * phoneme lookahead window, plus the current frame, plus one for
     * good measure? (FIXME: I don't remember why) */
    s->n_fast_hist = cmd_ln_int32_r(s->config, "-pl_window") + 2;
    s->hist = ckd_calloc(s->n_fast_hist, sizeof(*s->hist));
    /* s->f will be a rotating pointer into s->hist. */
    s->f = s->hist;
    for (i = 0; i < s->n_fast_hist; ++i) {
        int j, k, m;
        /* Top-N codewords for every codebook and feature. */
        s->hist[i].topn = ckd_calloc_3d(s->g->n_mgau, s->g->n_feat,
                                        s->max_topn, sizeof(ptm_topn_t));
        /* Initialize them to sane (yet arbitrary) defaults. */
        for (j = 0; j < s->g->n_mgau; ++j) {
            for (k = 0; k < s->g->n_feat; ++k) {
                for (m = 0; m < s->max_topn; ++m) {
                    s->hist[i].topn[j][k][m].cw = m;
                    s->hist[i].topn[j][k][m].score = WORST_DIST;
                }
            }
        }
        /* Active codebook mapping (just codebook, not features,
           at least not yet) */
        s->hist[i].mgau_active = bitvec_alloc(s->g->n_mgau);
        /* Start with them all on, prune them later. */
        bitvec_set_all(s->hist[i].mgau_active, s->g->n_mgau);
    }

    ps = (ps_mgau_t *)s;
    ps->vt = &ptm_mgau_funcs;
    return ps;
error_out:
    ptm_mgau_free(ps_mgau_base(s));
    return NULL;
}

ps_mgau_t *
ptm_mgau_copy(ps_mgau_t *other)
{
    ptm_mgau_t *others = (ptm_mgau_t *)other;
    ptm_mgau_t *s;
    ps_mgau_t *ps;
    int i;

    s = ckd_calloc(1, sizeof(*s));
    s->config = others->config;
    s->lmath = logmath_retain(others->lmath);
    s->lmath_8b = logmath_retain(others->lmath_8b);
    s->g = gauden_retain(others->g);
    s->n_sen = others->n_sen;
    s->ds_ratio = others->ds_ratio;
    s->max_topn = others->max_topn;
    s->s = sendump_retain(others->s);

    /* Allocate fast-match history buffers.  We need enough for the
     * phoneme lookahead window, plus the current frame, plus one for
     * good measure? (FIXME: I don't remember why) */
    /* FIXME: actually won't need this in the future when mgau copies
     * are used everywhere. */
    s->n_fast_hist = cmd_ln_int32_r(s->config, "-pl_window") + 2;
    s->hist = ckd_calloc(s->n_fast_hist, sizeof(*s->hist));
    /* s->f will be a rotating pointer into s->hist. */
    s->f = s->hist;
    for (i = 0; i < s->n_fast_hist; ++i) {
        int j, k, m;
        /* Top-N codewords for every codebook and feature. */
        s->hist[i].topn = ckd_calloc_3d(s->g->n_mgau, s->g->n_feat,
                                        s->max_topn, sizeof(ptm_topn_t));
        /* Initialize them to sane (yet arbitrary) defaults. */
        for (j = 0; j < s->g->n_mgau; ++j) {
            for (k = 0; k < s->g->n_feat; ++k) {
                for (m = 0; m < s->max_topn; ++m) {
                    s->hist[i].topn[j][k][m].cw = m;
                    s->hist[i].topn[j][k][m].score = WORST_DIST;
                }
            }
        }
        /* Active codebook mapping (just codebook, not features,
           at least not yet) */
        s->hist[i].mgau_active = bitvec_alloc(s->g->n_mgau);
        /* Start with them all on, prune them later. */
        bitvec_set_all(s->hist[i].mgau_active, s->g->n_mgau);
    }

    ps = (ps_mgau_t *)s;
    ps->vt = &ptm_mgau_funcs;
    return ps;
}

void
ptm_mgau_free(ps_mgau_t *ps)
{
    ptm_mgau_t *s = (ptm_mgau_t *)ps;

    logmath_free(s->lmath);
    logmath_free(s->lmath_8b);
    gauden_free(s->g);
    sendump_free(s->s);
    ckd_free(s);
}
