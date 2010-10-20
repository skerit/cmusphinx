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
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/sync_array.h>
#include <sphinxbase/sbthread.h>
#include <sphinxbase/ckd_alloc.h>

/* Local headers. */
#include "featbuf.h"

struct featbuf_s {
    sync_array_t *sa;
    cmd_ln_t *config;
    fe_t *fe;
    feat_t *fcb;
    mfcc_t *cepbuf;
    mfcc_t ***featbuf;
    int beginutt, endutt;
};

static int
featbuf_init_feat(featbuf_t *fb)
{
    fb->fcb = 
        feat_init(cmd_ln_str_r(fb->config, "-feat"),
                  cmn_type_from_str(cmd_ln_str_r(fb->config,"-cmn")),
                  cmd_ln_boolean_r(fb->config, "-varnorm"),
                  agc_type_from_str(cmd_ln_str_r(fb->config, "-agc")),
                  1, cmd_ln_int32_r(fb->config, "-ceplen"));
    if (fb->fcb == NULL)
        return -1;

    if (cmd_ln_str_r(fb->config, "-lda")) {
        E_INFO("Reading linear feature transformation from %s\n",
               cmd_ln_str_r(fb->config, "-lda"));
        if (feat_read_lda(fb->fcb,
                          cmd_ln_str_r(fb->config, "-lda"),
                          cmd_ln_int32_r(fb->config, "-ldadim")) < 0)
            return -1;
    }

    if (cmd_ln_str_r(fb->config, "-svspec")) {
        int32 **subvecs;
        E_INFO("Using subvector specification %s\n", 
               cmd_ln_str_r(fb->config, "-svspec"));
        if ((subvecs = parse_subvecs
             (cmd_ln_str_r(fb->config, "-svspec"))) == NULL)
            return -1;
        if ((feat_set_subvecs(fb->fcb, subvecs)) < 0)
            return -1;
    }

    if (cmd_ln_exists_r(fb->config, "-agcthresh")
        && 0 != strcmp(cmd_ln_str_r(fb->config, "-agc"), "none")) {
        agc_set_threshold(fb->fcb->agc_struct,
                          cmd_ln_float32_r(fb->config, "-agcthresh"));
    }

    if (fb->fcb->cmn_struct
        && cmd_ln_exists_r(fb->config, "-cmninit")) {
        char *c, *cc, *vallist;
        int32 nvals;

        vallist = ckd_salloc(cmd_ln_str_r(fb->config, "-cmninit"));
        c = vallist;
        nvals = 0;
        while (nvals < fb->fcb->cmn_struct->veclen
               && (cc = strchr(c, ',')) != NULL) {
            *cc = '\0';
            fb->fcb->cmn_struct->cmn_mean[nvals] = FLOAT2MFCC(atof(c));
            c = cc + 1;
            ++nvals;
        }
        if (nvals < fb->fcb->cmn_struct->veclen && *c != '\0') {
            fb->fcb->cmn_struct->cmn_mean[nvals] = FLOAT2MFCC(atof(c));
        }
        ckd_free(vallist);
    }

    return 0;
}

featbuf_t *
featbuf_init(cmd_ln_t *config)
{
    featbuf_t *fb;

    fb = ckd_calloc(1, sizeof(*fb));
    fb->config = cmd_ln_retain(config);
    fb->fe = fe_init_auto_r(config);
    if (featbuf_init_feat(fb) < 0)
        goto error_out;
    /* Only allocate one frame of features, not clear if it will be
     * faster to do more or not (we can only queue one frame at a time
     * for the moment) */
    fb->cepbuf = ckd_calloc(fe_get_output_size(fb->fe),
                            sizeof(*fb->cepbuf));
    fb->featbuf = feat_array_alloc(fb->fcb,
                                   feat_window_size(fb->fcb) + 1);
    return fb;
error_out:
    featbuf_free(fb);
    return NULL;
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
    ckd_free(fb->cepbuf);
    feat_array_free(fb->featbuf);
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

int
featbuf_wait(featbuf_t *fb, int fidx, int timeout, mfcc_t *out_frame)
{
    int s = timeout == -1 ? -1 : 0;
    int rc;

    /* Wait for the frame in sync array. */
    if ((rc = sync_array_wait(fb->sa, fidx, s, timeout)) < 0)
        return rc;

    /* Copy it. */
    sync_array_get(fb->sa, fidx, out_frame);

    return 0;
}

int
featbuf_release(featbuf_t *fb, int sidx, int eidx)
{
    return sync_array_release(fb->sa, sidx, eidx);
}

int
featbuf_start_utt(featbuf_t *fb)
{
    /* Reset the sync array. */
    sync_array_reset(fb->sa);

    /* Set utterance processing state. */
    fb->beginutt = TRUE;
    fb->endutt = FALSE;

    return -1;
}

int
featbuf_end_utt(featbuf_t *fb, int timeout)
{
    int nfr;

    /* Set utterance processing state. */
    fb->endutt = TRUE;

    /* Drain remaining frames from fb->fe. */
    if (fe_end_utt(fb->fe, fb->cepbuf, &nfr) < 0)
        return -1;

    /* Drain remaining frames from fb->fcb.*/
    if (featbuf_process_cep(fb, &fb->cepbuf, nfr, FALSE) < 0)
        return -1;

    /* Finalize. */
    sync_array_finalize(fb->sa);

    /* Wait for everybody to be done (FIXME: not quite sure how to do
     * that yet) */
    return 0;
}

int
featbuf_abort_utt(featbuf_t *fb)
{
    sync_array_force_quit(fb->sa);

    /* Wait for everybody to quit. */
    return 0;
}

int
featbuf_process_raw(featbuf_t *fb,
                    int16 const *raw,
                    size_t n_samps,
                    int full_utt)
{
    int16 const *rptr;

    /* Full utt, process it all (CMN, etc)... */
    if (full_utt) {
    }

    /* Write audio to log file. */

    /* Do fe_process_frames into our internal MFCC buffer until no
     * data remains. */
    rptr = raw;
    while (n_samps > 0) {
        int nframes = 1;
        if (fe_process_frames(fb->fe, &rptr, &n_samps,
                              &fb->cepbuf, &nframes) < 0)
            return -1;
        if (nframes)
            featbuf_process_cep(fb, &fb->cepbuf, 1, FALSE);
    }

    return 0;
}

int
featbuf_process_cep(featbuf_t *fb,
                    mfcc_t **cep,
                    size_t n_frames,
                    int full_utt)
{
    mfcc_t **cptr;

    /* Full utt, process it all (CMN, etc)... */
    if (full_utt) {
    }

    /* Write frames to log file. */

    /* Do s2mfc2feat_live on one frame at a time, appending the
     * resulting feature blocks. */
    cptr = cep;
    while (n_frames > 0) {
        int nfeat, i;
        int ncep = n_frames;
        nfeat = feat_s2mfc2feat_live(fb->fcb, cptr,
                                     &ncep, fb->beginutt, fb->endutt,
                                     fb->featbuf);
        if (fb->beginutt)
            fb->beginutt = FALSE;
        if (fb->endutt)
            assert(nfeat == n_frames);
        for (i = 0; i < nfeat; ++i) {
            if (featbuf_process_feat(fb, fb->featbuf[i]) < 0)
                return -1;
        }
        cptr += ncep;
        n_frames -= ncep;
    }

    return 0;
}

int
featbuf_process_feat(featbuf_t *fb,
                     mfcc_t **feat)
{
    /* This one is easy, at least... */
    return sync_array_append(fb->sa, feat[0]);
}
