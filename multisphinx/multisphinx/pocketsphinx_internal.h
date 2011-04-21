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

/* Other headers. */
#include <multisphinx/search.h>
#include <multisphinx/featbuf.h>
#include <multisphinx/acmod.h>
#include <multisphinx/dict.h>
#include <multisphinx/dict2pid.h>
#include <multisphinx/pocketsphinx.h>

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
    search_t *fwdtree;
    search_t *fwdflat;
};

#endif /* __POCKETSPHINX_INTERNAL_H__ */
