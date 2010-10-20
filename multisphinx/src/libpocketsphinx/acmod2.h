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
 * @file acmod2.h Acoustic model structures for PocketSphinx.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __ACMOD2_H__
#define __ACMOD2_H__

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/fe.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/bitvec.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "bin_mdef.h"
#include "tmat.h"
#include "hmm.h"

/**
 * Dummy senone score value for unintentionally active states.
 */
#define SENSCR_DUMMY 0x7fff

/**
 * Acoustic model parameter structure. 
 */
typedef struct ps_mgau_s ps_mgau_t;

typedef struct ps_mgaufuncs_s {
    char const *name;

    int (*frame_eval)(ps_mgau_t *mgau,
                      int16 *senscr,
                      uint8 *senone_active,
                      int32 n_senone_active,
                      mfcc_t ** feat,
                      int32 frame,
                      int32 compallsen);
    ps_mgau_t *(*copy)(ps_mgau_t *mgau);
    void (*free)(ps_mgau_t *mgau);
} ps_mgaufuncs_t;    

struct ps_mgau_s {
    ps_mgaufuncs_t *vt;  /**< vtable of mgau functions. */
    int frame_idx;       /**< frame counter. */
};

#define ps_mgau_base(mg) ((ps_mgau_t *)(mg))
#define ps_mgau_frame_eval(mg,senscr,senone_active,n_senone_active,feat,frame,compallsen) \
    (*ps_mgau_base(mg)->vt->frame_eval)                                 \
    (mg, senscr, senone_active, n_senone_active, feat, frame, compallsen)
#define ps_mgau_transform(mg, mllr)                                  \
    (*ps_mgau_base(mg)->vt->transform)(mg, mllr)
#define ps_mgau_free(mg)                                  \
    (*ps_mgau_base(mg)->vt->free)(mg)
#define ps_mgau_copy(mg)                                  \
    (*ps_mgau_base(mg)->vt->copy)(mg)

/**
 * Acoustic model structure.
 *
 * This object corresponds to the acoustic scoring portion of acoustic
 * evaluation.  Feature computation and buffering is done by
 * featbuf_t.  Since the two objects must agree on the parameters for
 * feature computation it is suggested that the same configuration
 * (cmd_ln_t) be used to initialize them.  Compatibility will be
 * checked at initialization time, however.
 */
struct acmod_s {
    int refcount;

    /* Global objects. */
    cmd_ln_t *config;          /**< Configuration. */
    logmath_t *lmath;          /**< Log-math computation. */
    glist_t strings;           /**< Temporary acoustic model filenames. */

    /* Model parameters: */
    bin_mdef_t *mdef;          /**< Model definition. */
    tmat_t *tmat;              /**< Transition matrices. */
    ps_mgau_t *mgau;           /**< Model parameters. */

    /* Senone scoring: */
    int16 *senone_scores;      /**< GMM scores for current frame. */
    bitvec_t *senone_active_vec; /**< Active GMMs in current frame. */
    uint8 *senone_active;      /**< Array of deltas to active GMMs. */
    int senscr_frame;          /**< Frame index for senone_scores. */
    int n_senone_active;       /**< Number of active GMMs. */
    int log_zero;              /**< Zero log-probability value. */

    /* Flags and counters: */
    int16 compallsen;   /**< Compute all senones? */
    int16 output_frame; /**< Index of next frame of dynamic features. */
};
typedef struct acmod_s acmod_t;

/**
 * Initialize an acoustic model.
 *
 * @param config A command-line object containing parameters.
 * @param lmath Global log-math parameters.
 * @param featbuf_t Feature buffer to obtain features from.
 * @return a newly initialized acmod_t, or NULL on failure.
 */
acmod_t *acmod_init(cmd_ln_t *config, logmath_t *lmath, featbuf_t *fb);

/**
 * Create a partial copy of an acoustic model.
 *
 * The new acoustic model shares the same parameters, but keeps its
 * own state, including acoustic feature buffers and the current frame
 * of evaluation.
 */
acmod_t *acmod_copy(acmod_t *acmod);

/**
 * Retain a pointer to an acoustic model.
 *
 * Unlike in acmod_copy() this is simply a reference to the same acoustic model.
 */
acmod_t *acmod_retain(acmod_t *acmod);


/**
 * Release a pointer to an acoustic model.
 */
int acmod_free(acmod_t *acmod);

/**
 * Score one frame of data.
 *
 * @param frame_idx Frame index to score.
 * @return Array of senone scores for this frame, or NULL on error.
 *         The data pointed to persists only until the next call to
 *         acmod_score() or acmod_advance().
 */
int16 const *acmod_score(acmod_t *acmod, int frame_idx);

/**
 * Get best score and senone index for current frame.
 */
int acmod_best_score(acmod_t *acmod, int *out_best_senid);

/**
 * Clear set of active senones.
 */
void acmod_clear_active(acmod_t *acmod);

/**
 * Activate senones associated with an HMM.
 */
void acmod_activate_hmm(acmod_t *acmod, hmm_t *hmm);

/**
 * Activate a single senone.
 */
#define acmod_activate_sen(acmod, sen) bitvec_set((acmod)->senone_active_vec, sen)

/**
 * Build active list from 
 */
int32 acmod_flags2list(acmod_t *acmod);

#endif /* __ACMOD2_H__ */
