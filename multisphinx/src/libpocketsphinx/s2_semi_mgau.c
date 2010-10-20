/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
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
#include "s2_semi_mgau.h"
#include "tied_mgau_common.h"
#include "posixwin32.h"

static ps_mgaufuncs_t s2_semi_mgau_funcs = {
    "s2_semi",
    &s2_semi_mgau_frame_eval,      /* frame_eval */
    &s2_semi_mgau_copy,            /* copy */
    &s2_semi_mgau_free             /* free */
};

struct vqFeature_s {
    int32 score; /* score or distance */
    int32 codeword; /* codeword (vector index) */
};

static void
eval_topn(s2_semi_mgau_t *s, int32 feat, mfcc_t *z)
{
    int i, ceplen;
    vqFeature_t *topn;

    topn = s->f[feat];
    ceplen = s->veclen[feat];

    for (i = 0; i < s->max_topn; i++) {
        mfcc_t *mean, diff, sqdiff, compl; /* diff, diff^2, component likelihood */
        vqFeature_t vtmp;
        mfcc_t *var, d;
        mfcc_t *obs;
        int32 cw, j;

        cw = topn[i].codeword;
        mean = s->means[feat][0] + cw * ceplen;
        var = s->vars[feat][0] + cw * ceplen;
        d = s->dets[feat][cw];
        obs = z;
        for (j = 0; j < ceplen; j++) {
            diff = *obs++ - *mean++;
            sqdiff = MFCCMUL(diff, diff);
            compl = MFCCMUL(sqdiff, *var);
            d = GMMSUB(d, compl);
            ++var;
        }
        topn[i].score = (int32)d;
        if (i == 0)
            continue;
        vtmp = topn[i];
        for (j = i - 1; j >= 0 && (int32)d > topn[j].score; j--) {
            topn[j + 1] = topn[j];
        }
        topn[j + 1] = vtmp;
    }
}

static void
eval_cb(s2_semi_mgau_t *s, int32 feat, mfcc_t *z)
{
    vqFeature_t *worst, *best, *topn;
    mfcc_t *mean;
    mfcc_t *var, *det, *detP, *detE;
    int32 i, ceplen;

    best = topn = s->f[feat];
    worst = topn + (s->max_topn - 1);
    mean = s->means[feat][0];
    var = s->vars[feat][0];
    det = s->dets[feat];
    detE = det + s->n_density;
    ceplen = s->veclen[feat];

    for (detP = det; detP < detE; ++detP) {
        mfcc_t diff, sqdiff, compl; /* diff, diff^2, component likelihood */
        mfcc_t d;
        mfcc_t *obs;
        vqFeature_t *cur;
        int32 cw, j;

        d = *detP;
        obs = z;
        cw = detP - det;
        for (j = 0; (j < ceplen) && (d >= worst->score); ++j) {
            diff = *obs++ - *mean++;
            sqdiff = MFCCMUL(diff, diff);
            compl = MFCCMUL(sqdiff, *var);
            d = GMMSUB(d, compl);
            ++var;
        }
        if (j < ceplen) {
            /* terminated early, so not in topn */
            mean += (ceplen - j);
            var += (ceplen - j);
            continue;
        }
        if ((int32)d < worst->score)
            continue;
        for (i = 0; i < s->max_topn; i++) {
            /* already there, so don't need to insert */
            if (topn[i].codeword == cw)
                break;
        }
        if (i < s->max_topn)
            continue;       /* already there.  Don't insert */
        /* remaining code inserts codeword and dist in correct spot */
        for (cur = worst - 1; cur >= best && (int32)d >= cur->score; --cur)
            memcpy(cur + 1, cur, sizeof(vqFeature_t));
        ++cur;
        cur->codeword = cw;
        cur->score = (int32)d;
    }
}

static void
mgau_dist(s2_semi_mgau_t * s, int32 frame, int32 feat, mfcc_t * z)
{
    eval_topn(s, feat, z);

    /* If this frame is skipped, do nothing else. */
    if (frame % s->ds_ratio)
        return;

    /* Evaluate the rest of the codebook (or subset thereof). */
    eval_cb(s, feat, z);
}

static int
mgau_norm(s2_semi_mgau_t *s, int feat)
{
    int32 norm;
    int j;

    /* Compute quantized normalizing constant. */
    norm = s->f[feat][0].score >> SENSCR_SHIFT;

    /* Normalize the scores, negate them, and clamp their dynamic range. */
    for (j = 0; j < s->max_topn; ++j) {
        s->f[feat][j].score = -((s->f[feat][j].score >> SENSCR_SHIFT) - norm);
        if (s->f[feat][j].score > MAX_NEG_ASCR)
            s->f[feat][j].score = MAX_NEG_ASCR;
        if (s->topn_beam[feat] && s->f[feat][j].score > s->topn_beam[feat])
            break;
    }
    return j;
}

static int32
get_scores_8b_feat_6(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3, *pid_cw4, *pid_cw5;

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->s->mixw[i][s->f[i][3].codeword];
    pid_cw4 = s->s->mixw[i][s->f[i][4].codeword];
    pid_cw5 = s->s->mixw[i][s->f[i][5].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw2[sen] + s->f[i][2].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw3[sen] + s->f[i][3].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw4[sen] + s->f[i][4].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw5[sen] + s->f[i][5].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_5(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3, *pid_cw4;

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->s->mixw[i][s->f[i][3].codeword];
    pid_cw4 = s->s->mixw[i][s->f[i][4].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw2[sen] + s->f[i][2].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw3[sen] + s->f[i][3].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw4[sen] + s->f[i][4].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_4(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3;

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->s->mixw[i][s->f[i][3].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw2[sen] + s->f[i][2].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw3[sen] + s->f[i][3].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_3(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2;

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->s->mixw[i][s->f[i][2].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);
        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw2[sen] + s->f[i][2].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_2(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1;

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;

        tmp = fast_logmath_add(s->lmath_8b, tmp,
                               pid_cw1[sen] + s->f[i][1].score);

        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_1(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0;

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        int32 tmp = pid_cw0[sen] + s->f[i][0].score;
        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat_any(s2_semi_mgau_t * s, int i, int topn,
                       int16 *senone_scores, uint8 *senone_active,
                       int32 n_senone_active)
{
    int32 j, k, l;

    for (l = j = 0; j < n_senone_active; j++) {
        int sen = senone_active[j] + l;
        uint8 *pid_cw;
        int32 tmp;
        pid_cw = s->s->mixw[i][s->f[i][0].codeword];
        tmp = pid_cw[sen] + s->f[i][0].score;
        for (k = 1; k < topn; ++k) {
            pid_cw = s->s->mixw[i][s->f[i][k].codeword];
            tmp = fast_logmath_add(s->lmath_8b, tmp,
                                   pid_cw[sen] + s->f[i][k].score);
        }
        senone_scores[sen] += tmp;
        l = sen;
    }
    return 0;
}

static int32
get_scores_8b_feat(s2_semi_mgau_t * s, int i, int topn,
                   int16 *senone_scores, uint8 *senone_active, int32 n_senone_active)
{
    switch (topn) {
    case 6:
        return get_scores_8b_feat_6(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 5:
        return get_scores_8b_feat_5(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 4:
        return get_scores_8b_feat_4(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 3:
        return get_scores_8b_feat_3(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 2:
        return get_scores_8b_feat_2(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 1:
        return get_scores_8b_feat_1(s, i, senone_scores,
                                    senone_active, n_senone_active);
    default:
        return get_scores_8b_feat_any(s, i, topn, senone_scores,
                                      senone_active, n_senone_active);
    }
}

static int32
get_scores_8b_feat_all(s2_semi_mgau_t * s, int i, int topn, int16 *senone_scores)
{
    int32 j, k;

    for (j = 0; j < s->n_sen; j++) {
        uint8 *pid_cw;
        int32 tmp;
        pid_cw = s->s->mixw[i][s->f[i][0].codeword];
        tmp = pid_cw[j] + s->f[i][0].score;
        for (k = 1; k < topn; ++k) {
            pid_cw = s->s->mixw[i][s->f[i][k].codeword];
            tmp = fast_logmath_add(s->lmath_8b, tmp,
                                   pid_cw[j] + s->f[i][k].score);
        }
        senone_scores[j] += tmp;
    }
    return 0;
}

static int32
get_scores_4b_feat_6(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3, *pid_cw4, *pid_cw5;
    uint8 w_den[6][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->s->mixw_cb[j] + s->f[i][1].score;
        w_den[2][j] = s->s->mixw_cb[j] + s->f[i][2].score;
        w_den[3][j] = s->s->mixw_cb[j] + s->f[i][3].score;
        w_den[4][j] = s->s->mixw_cb[j] + s->f[i][4].score;
        w_den[5][j] = s->s->mixw_cb[j] + s->f[i][5].score;
    }

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->s->mixw[i][s->f[i][3].codeword];
    pid_cw4 = s->s->mixw[i][s->f[i][4].codeword];
    pid_cw5 = s->s->mixw[i][s->f[i][5].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
            cw = pid_cw4[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[4][cw]);
            cw = pid_cw5[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[5][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
            cw = pid_cw4[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[4][cw]);
            cw = pid_cw5[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[5][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_5(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3, *pid_cw4;
    uint8 w_den[5][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->s->mixw_cb[j] + s->f[i][1].score;
        w_den[2][j] = s->s->mixw_cb[j] + s->f[i][2].score;
        w_den[3][j] = s->s->mixw_cb[j] + s->f[i][3].score;
        w_den[4][j] = s->s->mixw_cb[j] + s->f[i][4].score;
    }

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->s->mixw[i][s->f[i][3].codeword];
    pid_cw4 = s->s->mixw[i][s->f[i][4].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
            cw = pid_cw4[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[4][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
            cw = pid_cw4[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[4][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_4(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2, *pid_cw3;
    uint8 w_den[4][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->s->mixw_cb[j] + s->f[i][1].score;
        w_den[2][j] = s->s->mixw_cb[j] + s->f[i][2].score;
        w_den[3][j] = s->s->mixw_cb[j] + s->f[i][3].score;
    }

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->s->mixw[i][s->f[i][2].codeword];
    pid_cw3 = s->s->mixw[i][s->f[i][3].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
            cw = pid_cw3[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[3][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_3(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1, *pid_cw2;
    uint8 w_den[3][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->s->mixw_cb[j] + s->f[i][1].score;
        w_den[2][j] = s->s->mixw_cb[j] + s->f[i][2].score;
    }

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];
    pid_cw2 = s->s->mixw[i][s->f[i][2].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
            cw = pid_cw2[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[2][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_2(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0, *pid_cw1;
    uint8 w_den[2][16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[0][j] = s->s->mixw_cb[j] + s->f[i][0].score;
        w_den[1][j] = s->s->mixw_cb[j] + s->f[i][1].score;
    }

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];
    pid_cw1 = s->s->mixw[i][s->f[i][1].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] >> 4;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[0][cw];
            cw = pid_cw1[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp, w_den[1][cw]);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_1(s2_semi_mgau_t * s, int i,
                     int16 *senone_scores, uint8 *senone_active,
                     int32 n_senone_active)
{
    int32 j, l;
    uint8 *pid_cw0;
    uint8 w_den[16];

    /* Precompute scaled densities. */
    for (j = 0; j < 16; ++j) {
        w_den[j] = s->s->mixw_cb[j] + s->f[i][0].score;
    }

    pid_cw0 = s->s->mixw[i][s->f[i][0].codeword];

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;

        if (n & 1) {
            cw = pid_cw0[n/2] >> 4;
            tmp = w_den[cw];
        }
        else {
            cw = pid_cw0[n/2] & 0x0f;
            tmp = w_den[cw];
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat_any(s2_semi_mgau_t * s, int i, int topn,
                       int16 *senone_scores, uint8 *senone_active,
                       int32 n_senone_active)
{
    int32 j, k, l;

    for (l = j = 0; j < n_senone_active; j++) {
        int n = senone_active[j] + l;
        int tmp, cw;
        uint8 *pid_cw;
    
        pid_cw = s->s->mixw[i][s->f[i][0].codeword];
        if (n & 1)
            cw = pid_cw[n/2] >> 4;
        else
            cw = pid_cw[n/2] & 0x0f;
        tmp = s->s->mixw_cb[cw] + s->f[i][0].score;
        for (k = 1; k < topn; ++k) {
            pid_cw = s->s->mixw[i][s->f[i][k].codeword];
            if (n & 1)
                cw = pid_cw[n/2] >> 4;
            else
                cw = pid_cw[n/2] & 0x0f;
            tmp = fast_logmath_add(s->lmath_8b, tmp,
                                   s->s->mixw_cb[cw] + s->f[i][k].score);
        }
        senone_scores[n] += tmp;
        l = n;
    }
    return 0;
}

static int32
get_scores_4b_feat(s2_semi_mgau_t * s, int i, int topn,
                   int16 *senone_scores, uint8 *senone_active, int32 n_senone_active)
{
    switch (topn) {
    case 6:
        return get_scores_4b_feat_6(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 5:
        return get_scores_4b_feat_5(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 4:
        return get_scores_4b_feat_4(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 3:
        return get_scores_4b_feat_3(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 2:
        return get_scores_4b_feat_2(s, i, senone_scores,
                                    senone_active, n_senone_active);
    case 1:
        return get_scores_4b_feat_1(s, i, senone_scores,
                                    senone_active, n_senone_active);
    default:
        return get_scores_4b_feat_any(s, i, topn, senone_scores,
                                      senone_active, n_senone_active);
    }
}

static int32
get_scores_4b_feat_all(s2_semi_mgau_t * s, int i, int topn, int16 *senone_scores)
{
    int j, last_sen;

    j = 0;
    /* Number of senones is always even, but don't overrun if it isn't. */
    last_sen = s->n_sen & ~1;
    while (j < last_sen) {
        uint8 *pid_cw;
        int32 tmp0, tmp1;
        int k;

        pid_cw = s->s->mixw[i][s->f[i][0].codeword];
        tmp0 = s->s->mixw_cb[pid_cw[j/2] & 0x0f] + s->f[i][0].score;
        tmp1 = s->s->mixw_cb[pid_cw[j/2] >> 4] + s->f[i][0].score;
        for (k = 1; k < topn; ++k) {
            int32 w_den0, w_den1;

            pid_cw = s->s->mixw[i][s->f[i][k].codeword];
            w_den0 = s->s->mixw_cb[pid_cw[j/2] & 0x0f] + s->f[i][k].score;
            w_den1 = s->s->mixw_cb[pid_cw[j/2] >> 4] + s->f[i][k].score;
            tmp0 = fast_logmath_add(s->lmath_8b, tmp0, w_den0);
            tmp1 = fast_logmath_add(s->lmath_8b, tmp1, w_den1);
        }
        senone_scores[j++] += tmp0;
        senone_scores[j++] += tmp1;
    }
    return 0;
}

/*
 * Compute senone scores for the active senones.
 */
int32
s2_semi_mgau_frame_eval(ps_mgau_t *ps,
                        int16 *senone_scores,
                        uint8 *senone_active,
                        int32 n_senone_active,
			mfcc_t ** featbuf, int32 frame,
			int32 compallsen)
{
    s2_semi_mgau_t *s = (s2_semi_mgau_t *)ps;
    int i, topn_idx;

    memset(senone_scores, 0, s->n_sen * sizeof(*senone_scores));
    /* No bounds checking is done here, which just means you'll get
     * semi-random crap if you request a frame in the future or one
     * that's too far in the past. */
    topn_idx = frame % s->n_topn_hist;
    s->f = s->topn_hist[topn_idx];
    for (i = 0; i < s->n_feat; ++i) {
        /* For past frames this will already be computed. */
        if (frame >= ps_mgau_base(ps)->frame_idx) {
            vqFeature_t **lastf;
            if (topn_idx == 0)
                lastf = s->topn_hist[s->n_topn_hist-1];
            else
                lastf = s->topn_hist[topn_idx-1];
            memcpy(s->f[i], lastf[i], sizeof(vqFeature_t) * s->max_topn);
            mgau_dist(s, frame, i, featbuf[i]);
            s->topn_hist_n[topn_idx][i] = mgau_norm(s, i);
        }
        if (s->s->mixw_cb) {
            if (compallsen)
                get_scores_4b_feat_all(s, i, s->topn_hist_n[topn_idx][i], senone_scores);
            else
                get_scores_4b_feat(s, i, s->topn_hist_n[topn_idx][i], senone_scores,
                                   senone_active, n_senone_active);
        }
        else {
            if (compallsen)
                get_scores_8b_feat_all(s, i, s->topn_hist_n[topn_idx][i], senone_scores);
            else
                get_scores_8b_feat(s, i, s->topn_hist_n[topn_idx][i], senone_scores,
                                   senone_active, n_senone_active);
        }
    }

    return 0;
}

static int
split_topn(char const *str, uint8 *out, int nfeat)
{
    char *topn_list = ckd_salloc(str);
    char *c, *cc;
    int i, maxn;

    c = topn_list;
    i = 0;
    maxn = 0;
    while (i < nfeat && (cc = strchr(c, ',')) != NULL) {
        *cc = '\0';
        out[i] = atoi(c);
        if (out[i] > maxn) maxn = out[i];
        c = cc + 1;
        ++i;
    }
    if (i < nfeat && *c != '\0') {
        out[i] = atoi(c);
        if (out[i] > maxn) maxn = out[i];
        ++i;
    }
    while (i < nfeat)
        out[i++] = maxn;

    ckd_free(topn_list);
    return maxn;
}


ps_mgau_t *
s2_semi_mgau_init(acmod_t *acmod)
{
    s2_semi_mgau_t *s;
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
    /* Currently only a single codebook is supported. */
    if (s->g->n_mgau != 1)
        goto error_out;
    /* FIXME: maintaining pointers for convenience for now */
    s->means = s->g->mean[0];
    s->vars = s->g->var[0];
    s->dets = s->g->det[0];
    s->veclen = s->g->featlen;    
    /* Verify n_feat and veclen, against acmod. */
    s->n_feat = s->g->n_feat;
    if (s->n_feat != feat_dimension1(acmod->fcb)) {
        E_ERROR("Number of streams does not match: %d != %d\n",
                s->n_feat, feat_dimension(acmod->fcb));
        goto error_out;
    }
    for (i = 0; i < s->n_feat; ++i) {
        if (s->veclen[i] != feat_dimension2(acmod->fcb, i)) {
            E_ERROR("Dimension of stream %d does not match: %d != %d\n",
                    s->veclen[i], feat_dimension2(acmod->fcb, i));
            goto error_out;
        }
    }
    s->n_density = s->g->n_density;
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

    /* Determine top-N for each feature */
    s->topn_beam = ckd_calloc(s->n_feat, sizeof(*s->topn_beam));
    s->max_topn = cmd_ln_int32_r(s->config, "-topn");
    split_topn(cmd_ln_str_r(s->config, "-topn_beam"), s->topn_beam, s->n_feat);
    E_INFO("Maximum top-N: %d ", s->max_topn);
    E_INFOCONT("Top-N beams:");
    for (i = 0; i < s->n_feat; ++i) {
        E_INFOCONT(" %d", s->topn_beam[i]);
    }
    E_INFOCONT("\n");

    /* Top-N scores from recent frames */
    s->n_topn_hist = cmd_ln_int32_r(s->config, "-pl_window") + 2;
    s->topn_hist = (vqFeature_t ***)
        ckd_calloc_3d(s->n_topn_hist, s->n_feat, s->max_topn,
                      sizeof(***s->topn_hist));
    s->topn_hist_n = ckd_calloc_2d(s->n_topn_hist, s->n_feat,
                                   sizeof(**s->topn_hist_n));
    for (i = 0; i < s->n_topn_hist; ++i) {
        int j;
        for (j = 0; j < s->n_feat; ++j) {
            int k;
            for (k = 0; k < s->max_topn; ++k) {
                s->topn_hist[i][j][k].score = WORST_DIST;
                s->topn_hist[i][j][k].codeword = k;
            }
        }
    }

    ps = (ps_mgau_t *)s;
    ps->vt = &s2_semi_mgau_funcs;
    return ps;
error_out:
    s2_semi_mgau_free(ps_mgau_base(s));
    return NULL;
}

ps_mgau_t *
s2_semi_mgau_copy(ps_mgau_t *other)
{
    s2_semi_mgau_t *others = (s2_semi_mgau_t *)other;
    s2_semi_mgau_t *s;
    ps_mgau_t *ps;
    int i;

    s = ckd_calloc(1, sizeof(*s));
    s->config = others->config;
    s->lmath = logmath_retain(others->lmath);
    s->lmath_8b = logmath_retain(others->lmath_8b);
    s->g = gauden_retain(others->g);

    /* FIXME: maintaining pointers for convenience for now. */
    s->means = s->g->mean[0];
    s->vars = s->g->var[0];
    s->dets = s->g->det[0];
    s->veclen = s->g->featlen;    
    s->n_feat = s->g->n_feat;
    s->n_density = s->g->n_density;
    s->s = sendump_retain(others->s);
    s->n_sen = others->n_sen;
    s->ds_ratio = others->ds_ratio;
    s->max_topn = others->max_topn;

    /* Copy topn_beam as it can't be refcounted easily. */
    s->topn_beam = ckd_calloc(s->n_feat, sizeof(*s->topn_beam));
    memcpy(s->topn_beam, others->topn_beam,
           s->n_feat * sizeof(*s->topn_beam));

    /* Top-N scores from recent frames */
    /* FIXME: actually won't need this in the future when mgau copies
     * are used everywhere. */
    s->n_topn_hist = cmd_ln_int32_r(s->config, "-pl_window") + 2;
    s->topn_hist = (vqFeature_t ***)
        ckd_calloc_3d(s->n_topn_hist, s->n_feat, s->max_topn,
                      sizeof(***s->topn_hist));
    s->topn_hist_n = ckd_calloc_2d(s->n_topn_hist, s->n_feat,
                                   sizeof(**s->topn_hist_n));
    for (i = 0; i < s->n_topn_hist; ++i) {
        int j;
        for (j = 0; j < s->n_feat; ++j) {
            int k;
            for (k = 0; k < s->max_topn; ++k) {
                s->topn_hist[i][j][k].score = WORST_DIST;
                s->topn_hist[i][j][k].codeword = k;
            }
        }
    }

    ps = (ps_mgau_t *)s;
    ps->vt = &s2_semi_mgau_funcs;
    return ps;
}

void
s2_semi_mgau_free(ps_mgau_t *ps)
{
    s2_semi_mgau_t *s = (s2_semi_mgau_t *)ps;

    logmath_free(s->lmath);
    logmath_free(s->lmath_8b);
    gauden_free(s->g);
    sendump_free(s->s);
    ckd_free(s->topn_beam);
    ckd_free_2d(s->topn_hist_n);
    ckd_free_3d((void **)s->topn_hist);
    ckd_free(s);
}
