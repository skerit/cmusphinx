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
#include <sphinxbase/fe.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/bitvec.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "ps_mllr.h"
#include "bin_mdef.h"
#include "tmat.h"
#include "hmm.h"

/**
 * States in utterance processing.
 */
typedef enum acmod_state_e {
    ACMOD_IDLE,		/**< Not in an utterance. */
    ACMOD_STARTED,      /**< Utterance started, no data yet. */
    ACMOD_PROCESSING,   /**< Utterance in progress. */
    ACMOD_ENDED         /**< Utterance ended, still buffering. */
} acmod_state_t;

/**
 * Dummy senone score value for unintentionally active states.
 */
#define SENSCR_DUMMY 0x7fff

/**
 * Feature space linear transform structure.
 */
struct ps_mllr_s {
    int refcnt;     /**< Reference count. */
    int n_class;    /**< Number of MLLR classes. */
    int n_feat;     /**< Number of feature streams. */
    int *veclen;    /**< Length of input vectors for each stream. */
    float32 ****A;  /**< Rotation part of mean transformations. */
    float32 ***b;   /**< Bias part of mean transformations. */
    float32 ***h;   /**< Diagonal transformation of variances. */
    int32 *cb2mllr; /**< Mapping from codebooks to transformations. */
};

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
    int (*transform)(ps_mgau_t *mgau,
                     ps_mllr_t *mllr);
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
 * This object encapsulates all stages of acoustic processing, from
 * raw audio input to acoustic score output.  The reason for grouping
 * all of these modules together is that they all have to "agree" in
 * their parameterizations, and the configuration of the acoustic and
 * dynamic feature computation is completely dependent on the
 * parameters used to build the original acoustic model (which should
 * by now always be specified in a feat.params file).
 *
 * Because there is not a one-to-one correspondence from blocks of
 * input audio or frames of input features to frames of acoustic
 * scores (due to dynamic feature calculation), results may not be
 * immediately available after input, and the output results will not
 * correspond to the last piece of data input.
 */
struct acmod_s {
    int refcount;

    /* Global objects, not retained (FIXME: unsure if that's wise or not). */
    cmd_ln_t *config;          /**< Configuration. */
    logmath_t *lmath;          /**< Log-math computation. */
    glist_t strings;           /**< Temporary acoustic model filenames. */

    /* Feature computation: */
    fe_t *fe;                  /**< Acoustic feature computation. */
    feat_t *fcb;               /**< Dynamic feature computation. */

    /* Model parameters: */
    bin_mdef_t *mdef;          /**< Model definition. */
    tmat_t *tmat;              /**< Transition matrices. */
    ps_mgau_t *mgau;           /**< Model parameters. */
    ps_mllr_t *mllr;           /**< Speaker transformation. */

    /* Senone scoring: */
    int16 *senone_scores;      /**< GMM scores for current frame. */
    bitvec_t *senone_active_vec; /**< Active GMMs in current frame. */
    uint8 *senone_active;      /**< Array of deltas to active GMMs. */
    int senscr_frame;          /**< Frame index for senone_scores. */
    int n_senone_active;       /**< Number of active GMMs. */
    int log_zero;              /**< Zero log-probability value. */

    /* Utterance processing: */
    mfcc_t **mfc_buf;   /**< Temporary buffer of acoustic features. */
    mfcc_t ***feat_buf; /**< Temporary buffer of dynamic features. */
    FILE *rawfh;        /**< File for writing raw audio data. */
    FILE *mfcfh;        /**< File for writing acoustic feature data. */
    FILE *senfh;        /**< File for writing senone score data. */
    FILE *insenfh;	/**< Input senone score file. */
    long *framepos;     /**< File positions of recent frames in senone file. */

    /* A whole bunch of flags and counters: */
    uint8 state;        /**< State of utterance processing. */
    uint8 compallsen;   /**< Compute all senones? */
    uint8 grow_feat;    /**< Whether to grow feat_buf. */
    uint8 insen_swap;   /**< Whether to swap input senone score. */
    int16 output_frame; /**< Index of next frame of dynamic features. */
    int16 n_mfc_alloc;  /**< Number of frames allocated in mfc_buf */
    int16 n_mfc_frame;  /**< Number of frames active in mfc_buf */
    int16 mfc_outidx;   /**< Start of active frames in mfc_buf */
    int16 n_feat_alloc; /**< Number of frames allocated in feat_buf */
    int16 n_feat_frame; /**< Number of frames active in feat_buf */
    int16 feat_outidx;  /**< Start of active frames in feat_buf */
};
typedef struct acmod_s acmod_t;

/**
 * Initialize an acoustic model.
 *
 * @param config a command-line object containing parameters.  This
 *               pointer is not retained by this object.
 * @param lmath global log-math parameters.
 * @param fe a previously-initialized acoustic feature module to use,
 *           or NULL to create one automatically.  If this is supplied
 *           and its parameters do not match those in the acoustic
 *           model, this function will fail.  This pointer is not retained.
 * @param fcb a previously-initialized dynamic feature module to use,
 *           or NULL to create one automatically.  If this is supplied
 *           and its parameters do not match those in the acoustic
 *           model, this function will fail.  This pointer is not retained.
 * @return a newly initialized acmod_t, or NULL on failure.
 */
acmod_t *acmod_init(cmd_ln_t *config, logmath_t *lmath, fe_t *fe, feat_t *fcb);

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
 * Adapt acoustic model using a linear transform.
 *
 * FIXME: The semantics of this with respect to copied acoustic models
 * are not defined at the moment.
 *
 * @param mllr The new transform to use, or NULL to update the existing
 *              transform.  The decoder retains ownership of this pointer,
 *              so you should not attempt to free it manually.  Use
 *              ps_mllr_retain() if you wish to reuse it
 *              elsewhere.
 * @return The updated transform object for this decoder, or
 *         NULL on failure.
 */
ps_mllr_t *acmod_update_mllr(acmod_t *acmod, ps_mllr_t *mllr);

/**
 * Start logging senone scores to a filehandle.
 *
 * @param acmod Acoustic model object.
 * @param logfh Filehandle to log to.
 * @return 0 for success, <0 on error.
 */
int acmod_set_senfh(acmod_t *acmod, FILE *senfh);

/**
 * Start logging MFCCs to a filehandle.
 *
 * @param acmod Acoustic model object.
 * @param logfh Filehandle to log to.
 * @return 0 for success, <0 on error.
 */
int acmod_set_mfcfh(acmod_t *acmod, FILE *logfh);

/**
 * Start logging raw audio to a filehandle.
 *
 * @param acmod Acoustic model object.
 * @param logfh Filehandle to log to.
 * @return 0 for success, <0 on error.
 */
int acmod_set_rawfh(acmod_t *acmod, FILE *logfh);

/**
 * Mark the start of an utterance.
 */
int acmod_start_utt(acmod_t *acmod);

/**
 * Mark the end of an utterance.
 */
int acmod_end_utt(acmod_t *acmod);

/**
 * Get the number of frames of available feature data.
 *
 * @return Number of available frames of feature data.
 */
int acmod_available(acmod_t *acmod);

/**
 * Advance the frame index.
 *
 * This function moves to the next frame of input data.  Subsequent
 * calls to acmod_score() will return scores for that frame, until the
 * next call to acmod_advance().
 *
 * If no new frames of input data are available, it does nothing, and
 * returns -1.
 *
 * @return New frame index, or -1 if no frames are available.
 */
int acmod_advance(acmod_t *acmod);

/**
 * Get the current frame index.
 *
 * @return Index of current frame to be scored.
 */
int acmod_frame(acmod_t *acmod);

/**
 * Set up a senone score dump file for input.
 *
 * @param insenfh File handle of dump file
 * @return 0 for success, <0 for failure
 */
int acmod_set_insenfh(acmod_t *acmod, FILE *insenfh);

/**
 * Read one frame of scores from senone score dump file.
 *
 * @return Number of frames read or <0 on error.
 */
int acmod_read_scores(acmod_t *acmod);

/**
 * Score one frame of data.
 *
 * @param inout_frame_idx Input: frame index to score, or NULL
 *                        to obtain scores for the most recent frame.
 *                        Output: frame index corresponding to this
 *                        set of scores.
 * @return Array of senone scores for this frame, or NULL if no frame
 *         is available for scoring (such as if a frame index is
 *         requested that is not yet or no longer available).  The
 *         data pointed to persists only until the next call to
 *         acmod_score() or acmod_advance().
 */
int16 const *acmod_score(acmod_t *acmod,
                         int *inout_frame_idx);

/**
 * Write senone dump file header.
 */
int acmod_write_senfh_header(acmod_t *acmod, FILE *logfh);

/**
 * Write a frame of senone scores to a dump file.
 */
int acmod_write_scores(acmod_t *acmod, int n_active, uint8 const *active,
                       int16 const *senscr, FILE *senfh);


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

#endif /* __ACMOD_H__ */
