/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2010 Carnegie Mellon University.  All rights
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
 * @file featbuf.h Feature extraction and buffering for PocketSphinx.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

/* System headers. */

/* SphinxBase headers. */
#include <sphinxbase/sync_array.h>
#include <sphinxbase/sbthread.h>

/* Local headers. */
#include "featbuf.h"

struct featbuf_s {
    sync_array_t *sa;
    cmd_ln_t *config;
    fe_t *fe;
    feat_t *fcb;
};

featbuf_t *
featbuf_init(cmd_ln_t *config)
{
    featbuf_t *fb;

    fb = ckd_calloc(1, sizeof(*fb));

    return fb;
}

featbuf_t *
featbuf_retain(featbuf_t *fb)
{
    /* Piggyback on the refcount of the sync array. */
    sync_array_retain(fb->sa);
    return fb;
}

int
featbuf_free(featbuf_t *fb)
{
    int rc;

    if (fb == NULL)
        return 0;
    /* Piggyback on the refcount of the sync array. */
    if ((rc = sync_array_free(fb->sa)) > 0)
        return rc;
    cmd_ln_free_r(fb->config);
    fe_free(fb->fe);
    feat_free(fb->fcb);
    ckd_free(fb);
    return 0;
}

fe_t *
featbuf_get_fe(featbuf_t *fb)
{
    return fb->fe;
}

feat_t *
featbuf_get_fcb(featbuf_t *fb)
{
    return fb->fcb;
}

int
featbuf_next(featbuf_t *fb)
{
    return sync_array_next_idx(fb->sa);
}

mfcc_t **
featbuf_wait(featbuf_t *fb, int fidx, int timeout)
{
}

int
featbuf_release(featbuf_t *fb, int sidx, int eidx)
{
}

int
featbuf_start_utt(featbuf_t *fb)
{
}

int
featbuf_end_utt(featbuf_t *fb, int timeout)
{
}

int
featbuf_abort_utt(featbuf_t *fb)
{
}

int
featbuf_process_raw(featbuf_t *fb,
			int16 const *raw,
			size_t n_samps,
			int full_utt)
{
}

int
featbuf_process_cep(featbuf_t *fb,
			mfcc_t **cep,
			size_t n_frames,
			int full_utt)
{
}

int
featbuf_process_feat(featbuf_t *fb,
                     mfcc_t **feat)
{
}
