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
 * @file acmod.h Acoustic model structures for PocketSphinx.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __ACMOD_H__
#define __ACMOD_H__

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/glist.h>
#include <sphinxbase/fe.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/bitvec.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "featbuf.h"
#include "bin_mdef.h"

struct tmat_s;
struct hmm_s;

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
    featbuf_t *fb;             /**< Source of features. */
    feat_t *fcb;               /**< Feature parameters (belongs to @a fb) */

    /* Model parameters: */
    bin_mdef_t *mdef;          /**< Model definition. */
    struct tmat_s *tmat;              /**< Transition matrices. */
    ps_mgau_t *mgau;           /**< Model parameters. */

    /* Senone scoring: */
    mfcc_t ***feat_buf;        /**< Features for current frame. */
    int16 *senone_scores;      /**< GMM scores for current frame. */
    bitvec_t *senone_active_vec; /**< Active GMMs in current frame. */
    uint8 *senone_active;      /**< Array of deltas to active GMMs. */
    int n_senone_active;       /**< Number of active GMMs. */
    int log_zero;              /**< Zero log-probability value. */

    /* Flags and counters: */
    int output_frame;          /**< Index of next frame to score. */
    int compallsen;            /**< Compute all senone scores. */
    int eou;                   /**< At end of utterance input. */
    char *uttid;
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
 * Wait for a new frame of data to become available.
 *
 * This function is used by the first pass of recognition, which is
 * entirely driven by the acoustic model.
 *
 * If an error, timeout, or end-of-utterance condition is signalled
 * this function will return a negative number.
 *
 * @param acmod Acoustic model.
 * @param timeout Number of nanoseconds to wait, or -1 to wait indefinitely.
 * @return Index of next available frame, or <0 on timeout or end of
 *         utterance.
 */
int acmod_consumer_wait(acmod_t *acmod, int timeout);

/**
 * Get the number of frames processed so far.
 *
 * @return Next frame index, or number of frames processed.
 */
int acmod_frame(acmod_t *acmod);

/**
 * Calculate acoustic scores for a given frame index.
 *
 * This function does not block.  If the given frame index is not
 * available it will return NULL.
 *
 * The scores returned will only persist until the next call to
 * acmod_score().  It is not safe to call this function from multiple
 * threads.
 *
 * @param acmod Acoustic model.
 * @param frame_idx Index of requested frame.
 * @return Pointer to temporary senone scores.
 */
int16 const *acmod_score(acmod_t *acmod, int frame_idx);

/**
 * Relinquish interest in a frame index.
 *
 * Once the scores for a given frame index are no longer needed, that
 * frame should be released with this function in order to save memory.
 *
 * @param acmod Acoustic model.
 * @param frame_idx Index of frame to release.
 * @return 0, or <0 on error.
 */
int acmod_consumer_release(acmod_t *acmod, int frame_idx);

/**
 * Report an end-of-utterance condition.
 *
 * @return TRUE if an end-of-utterance has been signalled.
 */
int acmod_eou(acmod_t *acmod);

/**
 * Wait for a new utterance to start.
 *
 * This function resets the internal frame counter and other data
 * structures in the acmod_t, then blocks until a new utterance
 * starts.
 *
 * @param acmod Acoustic model.
 * @return 0, or <0 on timeout or error
 */
int acmod_consumer_start_utt(acmod_t *acmod, int timeout);

/**
 * Clean up after the end of an utterance.
 *
 * This function must be called after an end-of-utterance condition is
 * detected (usually by acmod_wait returning -1).  It releases all
 * remaining frames being waited on by the acoustic model.
 *
 * @return 0, or <0 on timeout or error
 */
int acmod_consumer_end_utt(acmod_t *acmod);

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
void acmod_activate_hmm(acmod_t *acmod, struct hmm_s *hmm);

/**
 * Activate a single senone.
 */
#define acmod_activate_sen(acmod, sen) bitvec_set((acmod)->senone_active_vec, sen)

/**
 * Build active list from 
 */
int32 acmod_flags2list(acmod_t *acmod);

#endif /* __ACMOD_H__ */
