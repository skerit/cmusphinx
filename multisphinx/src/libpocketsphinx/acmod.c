/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
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
 * @file acmod.c Acoustic model structures for PocketSphinx.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

/* System headers. */
#include <assert.h>
#include <string.h>

/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>
#include <sphinxbase/err.h>
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/byteorder.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/bio.h>

/* Local headers. */
#include "cmdln_macro.h"
#include "acmod.h"
#include "s2_semi_mgau.h"
#include "ptm_mgau.h"
#include "ms_mgau.h"
#include "hmm.h"
#include "tmat.h"

#ifndef WORDS_BIGENDIAN
#define WORDS_BIGENDIAN 1
#endif

static int
acmod_init_am(acmod_t *acmod)
{
    char const *mdeffn, *tmatfn;

    /* Read model definition. */
    if ((mdeffn = cmd_ln_str_r(acmod->config, "-mdef")) == NULL) {
        E_ERROR("Must specify -mdef or -hmm\n");
        return -1;
    }

    if ((acmod->mdef = bin_mdef_read(acmod->config, mdeffn)) == NULL) {
        E_ERROR("Failed to read model definition from %s\n", mdeffn);
        return -1;
    }

    /* Read transition matrices. */
    if ((tmatfn = cmd_ln_str_r(acmod->config, "-tmat")) == NULL) {
        E_ERROR("No tmat file specified\n");
        return -1;
    }
    acmod->tmat = tmat_init(tmatfn, acmod->lmath,
                            cmd_ln_float32_r(acmod->config, "-tmatfloor"),
                            TRUE);

    /* Read the acoustic models. */
    if ((cmd_ln_str_r(acmod->config, "-mean") == NULL)
        || (cmd_ln_str_r(acmod->config, "-var") == NULL)
        || (cmd_ln_str_r(acmod->config, "-tmat") == NULL)) {
        E_ERROR("No mean/var/tmat files specified\n");
        return -1;
    }

    if (cmd_ln_str_r(acmod->config, "-senmgau")) {
        E_INFO("Using general multi-stream GMM computation\n");
        acmod->mgau = ms_mgau_init(acmod->config, acmod->lmath, acmod->mdef);
        if (acmod->mgau == NULL)
            return -1;
    }
    else {
        E_INFO("Attempting to use SCHMM computation module\n");
        if ((acmod->mgau = s2_semi_mgau_init(acmod)) == NULL) {
            E_INFO("Attempting to use PTHMM computation module\n");
            if ((acmod->mgau = ptm_mgau_init(acmod)) == NULL) {
                E_INFO("Falling back to general multi-stream GMM computation\n");
                acmod->mgau = ms_mgau_init(acmod->config, acmod->lmath, acmod->mdef);
                if (acmod->mgau == NULL)
                    return -1;
            }
        }
    }

    return 0;
}

acmod_t *
acmod_init(cmd_ln_t *config, logmath_t *lmath, featbuf_t *fb)
{
    acmod_t *acmod;

    acmod = ckd_calloc(1, sizeof(*acmod));
    acmod->refcount = 1;
    acmod->config = cmd_ln_retain(config);
    acmod->lmath = logmath_retain(lmath);
    acmod->fb = featbuf_retain(fb);
    acmod->fcb = featbuf_get_fcb(acmod->fb);

    /* Load acoustic model parameters. */
    if (acmod_init_am(acmod) < 0)
        goto error_out;

    /* Senone computation stuff. */
    acmod->senone_scores = ckd_calloc(bin_mdef_n_sen(acmod->mdef),
                                                     sizeof(*acmod->senone_scores));
    acmod->senone_active_vec = bitvec_alloc(bin_mdef_n_sen(acmod->mdef));
    acmod->senone_active = ckd_calloc(bin_mdef_n_sen(acmod->mdef),
                                                     sizeof(*acmod->senone_active));
    acmod->log_zero = logmath_get_zero(acmod->lmath);
    acmod->compallsen = cmd_ln_boolean_r(config, "-compallsen");

    acmod->feat_buf = feat_array_alloc(acmod->fcb, 1);
    return acmod;

error_out:
    acmod_free(acmod);
    return NULL;
}

int
acmod_free(acmod_t *acmod)
{
    if (acmod == NULL)
        return 0;
    if (--acmod->refcount > 0)
        return acmod->refcount;

    ckd_free(acmod->senone_scores);
    ckd_free(acmod->senone_active_vec);
    ckd_free(acmod->senone_active);

    if (acmod->mdef)
        bin_mdef_free(acmod->mdef);
    if (acmod->tmat)
        tmat_free(acmod->tmat);
    if (acmod->mgau)
        ps_mgau_free(acmod->mgau);

    featbuf_free(acmod->fb);
    feat_array_free(acmod->feat_buf);
    logmath_free(acmod->lmath);
    cmd_ln_free_r(acmod->config);
    ckd_free(acmod);
    return 0;
}

acmod_t *
acmod_retain(acmod_t *acmod)
{
    ++acmod->refcount;
    return acmod;
}

acmod_t *
acmod_copy(acmod_t *other)
{
    acmod_t *acmod;

    acmod = ckd_calloc(1, sizeof(*acmod));
    acmod->refcount = 1;
    acmod->config = cmd_ln_retain(other->config);
    acmod->lmath = logmath_retain(other->lmath);
    acmod->mdef = bin_mdef_retain(other->mdef);
    acmod->tmat = tmat_retain(other->tmat);
    acmod->mgau = ps_mgau_copy(other->mgau);
    acmod->fb = featbuf_retain(other->fb);
    acmod->fcb = other->fcb; /* Implicitly retained with fb, I think */

    /* Senone computation stuff. */
    acmod->senone_scores = ckd_calloc(bin_mdef_n_sen(acmod->mdef),
                                                     sizeof(*acmod->senone_scores));
    acmod->senone_active_vec = bitvec_alloc(bin_mdef_n_sen(acmod->mdef));
    acmod->senone_active = ckd_calloc(bin_mdef_n_sen(acmod->mdef),
                                                     sizeof(*acmod->senone_active));
    acmod->log_zero = logmath_get_zero(acmod->lmath);
    acmod->compallsen = cmd_ln_boolean_r(acmod->config, "-compallsen");

    acmod->feat_buf = feat_array_alloc(acmod->fcb, 1);

    return acmod;
}

int
acmod_consumer_wait(acmod_t *acmod, int timeout)
{
    int rv;

    if ((rv = featbuf_consumer_wait(acmod->fb, acmod->output_frame,
                                    timeout, acmod->feat_buf[0][0])) < 0) {
        E_INFO("EOU in frame %d\n", acmod->output_frame);
        /* This means end of utterance. */
        acmod->eou = TRUE;
        return rv;
    }
    return acmod->output_frame++;
}

int16 const *
acmod_score(acmod_t *acmod, int frame_idx)
{
    /* Obtain the frame to be scored. */
    if (featbuf_consumer_wait(acmod->fb, frame_idx,
                              0, acmod->feat_buf[0][0]) < 0)
        return NULL;

    /* Build active senone list. */
    acmod_flags2list(acmod);

    /* Generate scores for the next available frame */
    ps_mgau_frame_eval(acmod->mgau,
                       acmod->senone_scores,
                       acmod->senone_active,
                       acmod->n_senone_active,
                       acmod->feat_buf[0],
                       frame_idx,
                       acmod->compallsen);

    return acmod->senone_scores;
}

int
acmod_consumer_release(acmod_t *acmod, int frame_idx)
{
    return featbuf_consumer_release(acmod->fb, frame_idx, frame_idx + 1);
}

int
acmod_eou(acmod_t *acmod)
{
    return acmod->eou;
}

int
acmod_frame(acmod_t *acmod)
{
    return acmod->output_frame;
}

int
acmod_consumer_start_utt(acmod_t *acmod, int timeout)
{
    int rc;

    if ((rc = featbuf_consumer_start_utt(acmod->fb, timeout)) < 0) {
        return rc;
    }
    
    E_INFO("Finished waiting for start of utt\n");
    acmod->output_frame = 0;
    acmod->eou = FALSE;
    acmod->uttid = featbuf_uttid(acmod->fb);

    return 0;
}

int
acmod_consumer_end_utt(acmod_t *acmod)
{
    featbuf_consumer_end_utt(acmod->fb, acmod->output_frame);
    acmod->eou = TRUE;

    return 0;
}

int
acmod_best_score(acmod_t *acmod, int *out_best_senid)
{
    int i, best;

    best = SENSCR_DUMMY;
    if (acmod->compallsen) {
        for (i = 0; i < bin_mdef_n_sen(acmod->mdef); ++i) {
            if (acmod->senone_scores[i] < best) {
                best = acmod->senone_scores[i];
                *out_best_senid = i;
            }
        }
    }
    else {
        int16 *senscr;
        senscr = acmod->senone_scores;
        for (i = 0; i < acmod->n_senone_active; ++i) {
            senscr += acmod->senone_active[i];
            if (*senscr < best) {
                best = *senscr;
                *out_best_senid = i;
            }
        }
    }
    return best;
}


void
acmod_clear_active(acmod_t *acmod)
{
    if (acmod->compallsen)
        return;
    bitvec_clear_all(acmod->senone_active_vec, bin_mdef_n_sen(acmod->mdef));
    acmod->n_senone_active = 0;
}

#define MPX_BITVEC_SET(a,h,i)                                   \
    if (hmm_mpx_ssid(h,i) != BAD_SSID)                          \
        bitvec_set((a)->senone_active_vec, hmm_mpx_senid(h,i))
#define NONMPX_BITVEC_SET(a,h,i)                                        \
    bitvec_set((a)->senone_active_vec,                                  \
               hmm_nonmpx_senid(h,i))

void
acmod_activate_hmm(acmod_t *acmod, hmm_t *hmm)
{
    int i;

    if (acmod->compallsen)
        return;
    if (hmm_is_mpx(hmm)) {
        switch (hmm_n_emit_state(hmm)) {
        case 5:
            MPX_BITVEC_SET(acmod, hmm, 4);
            MPX_BITVEC_SET(acmod, hmm, 3);
        case 3:
            MPX_BITVEC_SET(acmod, hmm, 2);
            MPX_BITVEC_SET(acmod, hmm, 1);
            MPX_BITVEC_SET(acmod, hmm, 0);
            break;
        default:
            for (i = 0; i < hmm_n_emit_state(hmm); ++i) {
                MPX_BITVEC_SET(acmod, hmm, i);
            }
        }
    }
    else {
        switch (hmm_n_emit_state(hmm)) {
        case 5:
            NONMPX_BITVEC_SET(acmod, hmm, 4);
            NONMPX_BITVEC_SET(acmod, hmm, 3);
        case 3:
            NONMPX_BITVEC_SET(acmod, hmm, 2);
            NONMPX_BITVEC_SET(acmod, hmm, 1);
            NONMPX_BITVEC_SET(acmod, hmm, 0);
            break;
        default:
            for (i = 0; i < hmm_n_emit_state(hmm); ++i) {
                NONMPX_BITVEC_SET(acmod, hmm, i);
            }
        }
    }
}

int32
acmod_flags2list(acmod_t *acmod)
{
    int32 w, l, n, b, total_dists, total_words, extra_bits;
    bitvec_t *flagptr;

    total_dists = bin_mdef_n_sen(acmod->mdef);
    if (acmod->compallsen) {
        acmod->n_senone_active = total_dists;
        return total_dists;
    }
    total_words = total_dists / BITVEC_BITS;
    extra_bits = total_dists % BITVEC_BITS;
    w = n = l = 0;
    for (flagptr = acmod->senone_active_vec; w < total_words; ++w, ++flagptr) {
        if (*flagptr == 0)
            continue;
        for (b = 0; b < BITVEC_BITS; ++b) {
            if (*flagptr & (1UL << b)) {
                int32 sen = w * BITVEC_BITS + b;
                int32 delta = sen - l;
                /* Handle excessive deltas "lossily" by adding a few
                   extra senones to bridge the gap. */
                while (delta > 255) {
                    acmod->senone_active[n++] = 255;
                    delta -= 255;
                }
                acmod->senone_active[n++] = delta;
                l = sen;
            }
        }
    }

    for (b = 0; b < extra_bits; ++b) {
        if (*flagptr & (1UL << b)) {
            int32 sen = w * BITVEC_BITS + b;
            int32 delta = sen - l;
            /* Handle excessive deltas "lossily" by adding a few
               extra senones to bridge the gap. */
            while (delta > 255) {
                acmod->senone_active[n++] = 255;
                delta -= 255;
            }
            acmod->senone_active[n++] = delta;
            l = sen;
        }
    }

    acmod->n_senone_active = n;
    E_DEBUG(1, ("acmod_flags2list: %d active in frame %d\n",
                acmod->n_senone_active, acmod->output_frame));
    return n;
}
