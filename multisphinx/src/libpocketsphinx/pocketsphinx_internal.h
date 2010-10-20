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
 * @file pocketsphinx_internal.h Internal implementation of
 * PocketSphinx decoder.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __POCKETSPHINX_INTERNAL_H__
#define __POCKETSPHINX_INTERNAL_H__

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/fe.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/profile.h>

/* Local headers. */
#include "pocketsphinx.h"
#include "ps_search.h"
#include "featbuf.h"
#include "acmod.h"
#include "dict.h"
#include "dict2pid.h"

/**
 * V-table for segmentation iterators.
 */
typedef struct ps_segfuncs_s {
    ps_seg_t *(*seg_next)(ps_seg_t *seg);
    void (*seg_free)(ps_seg_t *seg);
} ps_segfuncs_t;

/**
 * Base structure for hypothesis segmentation iterator.
 */
struct ps_seg_s {
    ps_segfuncs_t *vt;     /**< V-table of seg methods */
    ps_search_t *search;   /**< Search object from whence this came */
    char const *word;      /**< Word string (pointer into dictionary hash) */
    int16 sf;                /**< Start frame. */
    int16 ef;                /**< End frame. */
    int32 ascr;            /**< Acoustic score. */
    int32 lscr;            /**< Language model score. */
    int32 prob;            /**< Log posterior probability. */
    /* This doesn't need to be 32 bits, so once the scores above are
     * reduced to 16 bits (or less!), this will be too. */
    int32 lback;           /**< Language model backoff. */
    /* Not sure if this should be here at all. */
    float32 lwf;           /**< Language weight factor (for second-pass searches) */
};

#define ps_search_seg_next(seg) (*(seg->vt->seg_next))(seg)
#define ps_search_seg_free(s) (*(seg->vt->seg_free))(seg)


/**
 * Decoder object - implements the user-visible API.
 */
struct ps_decoder_s {
    int refcount;      /**< Reference count. */
    cmd_ln_t *config;  /**< Configuration. */

    /* Utterance-processing related stuff. */
    uint32 uttno;       /**< Utterance counter. */
    char *uttid;        /**< Utterance ID for current utterance. */
    ptmr_t perf;        /**< Performance counter for all of decoding. */
    uint32 n_frame;     /**< Total number of frames processed. */
    featbuf_t *fb;      /**< Feature buffer. */
    logmath_t *lmath;   /**< Global log-math computation. */
    acmod_t *acmod;     /**< Initial, global acoustic model. */
    
    /* Search modules (each of which has its own thread and its own
     * acmod_t, which may be cloned). */
    /* FIXME: Currently the fwdtree->fwdflat topology is hardwired,
     * this will change very soon. */
    ps_search_t *fwdtree;
    ps_search_t *fwdflat;

};

#endif /* __POCKETSPHINX_INTERNAL_H__ */
