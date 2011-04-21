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
#include <config.h>

/* SphinxBase headers. */
#include <sphinxbase/sync_array.h>
#include <sphinxbase/sbthread.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/byteorder.h>

/* Local headers. */
#include "featbuf.h"

/* We use normal C conditionals for this not preprocessor ones. */
#ifndef WORDS_BIGENDIAN
#define WORDS_BIGENDIAN 0
#endif

struct featbuf_s {
    int refcount;
    sync_array_t *sa;
    cmd_ln_t *config;
    fe_t *fe;
    feat_t *fcb;
    mfcc_t *cepbuf;
    mfcc_t ***featbuf;
    int beginutt, endutt;
    FILE *mfcfh;
    FILE *rawfh;

    sbsem_t *release;
    sbsem_t *start;

    /**
     * Flag signifying that featbuf_cancel() was called.  Reset by
     * featbuf_start_utt().  */
    int canceled;
    char *uttid;
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
    fb->refcount = 1;
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
    /* Each element is a complete (flattened) frame of features. */
    fb->sa = sync_array_init(0,
                             feat_dimension(fb->fcb)
                             * sizeof(mfcc_t));
    fb->start = sbsem_init("featbuf:start",0);
    fb->release = sbsem_init("featbuf:release",0);
    return fb;
error_out:
    featbuf_free(fb);
    return NULL;
}

featbuf_t *
featbuf_retain(featbuf_t *fb)
{
    ++fb->refcount;
    sync_array_retain(fb->sa);
    return fb;
}

int
featbuf_free(featbuf_t *fb)
{
    if (fb == NULL)
        return 0;
    sync_array_free(fb->sa);
    if (--fb->refcount > 0)
        return fb->refcount;

    /* Not refcounting these things internally. */
    fe_free(fb->fe);
    feat_free(fb->fcb);
    /* Not really sure why we can't do this. */
    /* cmd_ln_free_r(fb->config); */

    /* Non-refcounted things. */
    sbsem_free(fb->release);
    sbsem_free(fb->start);
    ckd_free(fb->cepbuf);
    feat_array_free(fb->featbuf);
    if (fb->mfcfh)
        fclose(fb->mfcfh);
    if (fb->rawfh)
        fclose(fb->rawfh);
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
featbuf_consumer_start_utt(featbuf_t *fb, int timeout)
{
    int s = (timeout == -1) ? -1 : 0;
    int rc;

    /* Wait for the semaphore to be zero then record this thread as
     * started. */
    if ((rc = sbsem_down(fb->start, s, timeout)) < 0)
        return rc;
    if (fb->canceled)
        return -1;
    return 0;
}

int
featbuf_consumer_wait(featbuf_t *fb, int fidx, int timeout, mfcc_t *out_frame)
{
    int s = timeout == -1 ? -1 : 0;
    int rc;

    /* Wait for the frame in sync array.  <0 means timeout or end of
     * utterance.*/
    if ((rc = sync_array_wait(fb->sa, fidx, s, timeout)) < 0)
        return rc;

    /* Copy it. */
    sync_array_get(fb->sa, fidx, out_frame);

    return 0;
}

int
featbuf_consumer_release(featbuf_t *fb, int sidx, int eidx)
{
    int rv;

    if (eidx == -1)
        eidx = sync_array_next_idx(fb->sa);
    if ((rv = sync_array_release(fb->sa, sidx, eidx)) < 0)
        return rv;
    return rv;
}

int
featbuf_consumer_end_utt(featbuf_t *fb, int sidx)
{
    int rv;

    if ((rv = featbuf_consumer_release(fb, sidx, -1)) < 0)
        return rv;
    /* Record this thread as having finished. */
    sbsem_up(fb->release);
    return rv;
}

int
featbuf_producer_start_utt(featbuf_t *fb, char *uttid)
{
    /* Reset the sync array. */
    sync_array_reset(fb->sa);

    /* Set utterance processing state. */
    fb->beginutt = TRUE;
    fb->endutt = FALSE;
    fb->uttid = uttid;

    fe_start_utt(fb->fe);

    /* Signal any consumers. */
    E_INFO("Setting canceled = FALSE\n");
    fb->canceled = FALSE;
    /* Allow refcount - 1 threads to start. */
    E_INFO("setting fb->start to %d\n", fb->refcount - 1);
    sbsem_set(fb->start, fb->refcount - 1);

    return 0;
}

int
featbuf_producer_end_utt(featbuf_t *fb)
{
    int nfr, i, rc, nth;
    size_t last_idx;

    /* Set utterance processing state. */
    fb->endutt = TRUE;

    /* Drain remaining frames from fb->fe. */
    if (fe_end_utt(fb->fe, fb->cepbuf, &nfr) < 0)
        return -1;

    /* Drain remaining frames from fb->fcb.*/
    if (featbuf_producer_process_cep(fb, &fb->cepbuf, nfr, FALSE) < 0)
        return -1;

    /* Close out log files. */
    if (fb->mfcfh) {
        int32 outlen, rv;
        outlen = (ftell(fb->mfcfh) - 4) / 4;
        if (!WORDS_BIGENDIAN)
            SWAP_INT32(&outlen);
        /* Try to seek and write */
        if ((rv = fseek(fb->mfcfh, 0, SEEK_SET)) == 0) {
            fwrite(&outlen, 4, 1, fb->mfcfh);
        }
        fclose(fb->mfcfh);
        fb->mfcfh = NULL;
    }
    if (fb->rawfh) {
        fclose(fb->rawfh);
        fb->rawfh = NULL;
    }

    /* Figure out how many threads we are going to be waiting for
     * before we finalize (since they may exit after that...) */
    nth = fb->refcount - 1;

    /* Finalize. */
    E_INFO("Finalizing frame array\n");
    last_idx = sync_array_finalize(fb->sa);

    /* Wait for everybody to be done. */
    for (i = 0; i < nth; ++i)
        if ((rc = sbsem_down(fb->release, -1, -1)) < 0)
            return rc;
    
    return 0;
}

int
featbuf_producer_shutdown(featbuf_t *fb)
{
    /* Wake up anybody waiting for an utterance, but first set a flag
     * that makes featbuf_wait_utt() fail. */
    E_INFO("Setting canceled = TRUE\n");
    fb->canceled = TRUE;
    E_INFO("setting fb->start to %d\n", fb->refcount - 1);
    sbsem_set(fb->start, fb->refcount - 1);
    return 0;
}

static int
featbuf_process_full_cep(featbuf_t *fb,
                         mfcc_t **cep,
                         int n_frames)
{
    mfcc_t ***featbuf;
    int nfr, i;

    /* Make a whole utterance of dynamic features. */
    featbuf = feat_array_alloc(fb->fcb, n_frames);
    nfr = feat_s2mfc2feat_live(fb->fcb, cep, &n_frames,
                               TRUE, TRUE, featbuf);

    /* Queue them. */
    for (i = 0; i < nfr; ++i) {
        if (featbuf_producer_process_feat(fb, featbuf[i]) < 0) {
            feat_array_free(featbuf);
            return -1;
        }
    }

    /* Don't worry about this extra allocation, the overhead of doing
     * this is almost certainly trivial compared to the cost of
     * decoding (plus CMN works better). */
    feat_array_free(featbuf);
    return nfr;
}

static int
featbuf_process_full_raw(featbuf_t *fb,
                         int16 const *raw,
                         size_t n_samps)
{
    mfcc_t **cepbuf = NULL;
    int nfr, ntail;

    if (fe_process_frames(fb->fe, NULL, &n_samps, NULL, &nfr) < 0)
        return -1;

    cepbuf = ckd_calloc_2d(nfr + 1, fe_get_output_size(fb->fe),
                           sizeof(**cepbuf));
    fe_start_utt(fb->fe);
    if (fe_process_frames(fb->fe, &raw, &n_samps,
                          cepbuf, &nfr) < 0)
        goto error_out;
    fe_end_utt(fb->fe, cepbuf[nfr], &ntail);
    nfr += ntail;

    if (featbuf_process_full_cep(fb, cepbuf, nfr) < 0)
        goto error_out;
    ckd_free_2d(cepbuf);
    return nfr;
error_out:
    ckd_free_2d(cepbuf);
    return -1;
}

int
featbuf_producer_process_raw(featbuf_t *fb,
                             int16 const *raw,
                             size_t n_samps,
                             int full_utt)
{
    int16 const *rptr;
    int total_nfr;

    /* Write audio to log file. */
    if (fb->rawfh)
        fwrite(raw, 2, n_samps, fb->rawfh);

    /* Full utt, process it all (CMN, etc)... */
    if (full_utt)
        return featbuf_process_full_raw(fb, raw, n_samps);

    /* Do fe_process_frames into our internal MFCC buffer until no
     * data remains. */
    rptr = raw;
    total_nfr = 0;
    while (n_samps > 0) {
        int nframes = 1;
        if (fe_process_frames(fb->fe, &rptr, &n_samps,
                              &fb->cepbuf, &nframes) < 0)
            return -1;
        if (nframes)
            featbuf_producer_process_cep(fb, &fb->cepbuf, 1, FALSE);
        total_nfr += nframes;
    }

    return total_nfr;
}

static int
featbuf_log_mfc(featbuf_t *fb,
                mfcc_t **cep, int n_frames)
{
    int i, n;
    int32 *ptr = (int32 *)cep[0];

    n = n_frames * feat_cepsize(fb->fcb);
    /* Swap bytes. */
    if (!WORDS_BIGENDIAN) {
        for (i = 0; i < (n * sizeof(mfcc_t)); ++i) {
            SWAP_INT32(ptr + i);
        }
    }
    /* Write features. */
    if (fwrite(cep[0], sizeof(mfcc_t), n, fb->mfcfh) != n) {
        E_ERROR_SYSTEM("Failed to write %d values to log file", n);
    }

    /* Swap them back. */
    if (!WORDS_BIGENDIAN) {
        for (i = 0; i < (n * sizeof(mfcc_t)); ++i) {
            SWAP_INT32(ptr + i);
        }
    }
    return 0;
}

int
featbuf_producer_process_cep(featbuf_t *fb,
                             mfcc_t **cep,
                             size_t n_frames,
                             int full_utt)
{
    mfcc_t **cptr;
    int out_nframes = 0;

    /* Write frames to log file. */
    if (fb->mfcfh)
        featbuf_log_mfc(fb, cep, n_frames);

    /* Full utt, process it all (CMN, etc)... */
    if (full_utt)
        return featbuf_process_full_cep(fb, cep, n_frames);

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
        for (i = 0; i < nfeat; ++i) {
            if (featbuf_producer_process_feat(fb, fb->featbuf[i]) < 0)
                return -1;
        }
        cptr += ncep;
        n_frames -= ncep;
        out_nframes += nfeat;
    }

    return out_nframes;
}

int
featbuf_producer_process_feat(featbuf_t *fb,
                              mfcc_t **feat)
{
    /* This one is easy, at least... */
    if (sync_array_append(fb->sa, feat[0]) < 0)
        return -1;
    return 1;
}

int
featbuf_set_mfcfh(featbuf_t *fb, FILE *logfh)
{
    int rv = 0;

    if (fb->mfcfh)
        fclose(fb->mfcfh);
    fb->mfcfh = logfh;
    fwrite(&rv, 4, 1, fb->mfcfh);
    return rv;
}

int
featbuf_set_rawfh(featbuf_t *fb, FILE *logfh)
{
    if (fb->rawfh)
        fclose(fb->rawfh);
    fb->rawfh = logfh;
    return 0;
}

char *
featbuf_uttid(featbuf_t *fb)
{
    return fb->uttid;
}

int
featbuf_get_window_start(featbuf_t *fb)
{
    return sync_array_available(fb->sa);
}

int
featbuf_get_window_end(featbuf_t *fb)
{
    return sync_array_next_idx(fb->sa);
}
