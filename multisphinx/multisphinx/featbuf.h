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

#ifndef __FEATBUF_H__
#define __FEATBUF_H__

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/fe.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/err.h>

/* Local headers. */

typedef struct featbuf_s featbuf_t;

/**
 * Create and initialize a new feature buffer.
 */
featbuf_t *featbuf_init(cmd_ln_t *config);

/**
 * Retain a pointer to a feature buffer.
 */
featbuf_t *featbuf_retain(featbuf_t *fb);

/**
 * Release a pointer to a feature buffer.
 */
int featbuf_free(featbuf_t *fb);

/**
 * Get a pointer to the feature extraction component.
 *
 * @return Feature extractor, retained by @a fb, DO NOT FREE.
 */
fe_t *featbuf_get_fe(featbuf_t *fb);

/**
 * Get a pointer to the feature computation component.
 *
 * @return Feature computer, retained by @a fb, DO NOT FREE.
 */
feat_t *featbuf_get_fcb(featbuf_t *fb);

/**
 * Wait for the beginning of an utterance.
 *
 * If an utterance is already in progress this returns immediately.
 *
 * @param fb Feature buffer.
 * @param timeout Maximum time to wait, in nanoseconds, or -1 to wait forever.
 * @return 0, or <0 on timeout or failure.
 */
int featbuf_consumer_start_utt(featbuf_t *fb, int timeout);

/**
 * Get the index of the next frame to become available.
 *
 * This is also the same as the number of frames processed so far.
 *
 * @param fb Feature buffer.
 * @return Index of the next frame to become available.
 */
int featbuf_next(featbuf_t *fb);

/**
 * Wait for a frame and its successors to become available.
 *
 * This function blocks for the requested timeout, or until the
 * requested frame becomes available.  If the timeout is reached it
 * will return -1.  It will also return -1 in the case where the
 * end of the utterance has been established by featbuf_end_utt() and
 * the frame requested is beyond the end of the utterance.
 *
 * In effect there is a semaphore attached to each frame index, and
 * this can be thought of as the "down" operation for a frame index.
 *
 * Calling this with a non-zero timeout on the same thread which calls
 * the data processing functions below is a Bad Idea, but that should
 * be obvious.
 *
 * For reasons of thread safety it is not possible to return a pointer
 * to the actual frame data.  To avoid memory allocation the caller
 * must allocate sufficient space for a single frame of dynamic
 * feature data.  This can be accomplished with feat_array_alloc().
 * The pointer to pass here is the first stream of the given frame,
 * e.g. if you have allocated mfcc ***x, pass x[frame][0].
 *
 * @param fb Feature buffer.
 * @param fidx Index of frame requested.
 * @param timeout Maximum time to wait, in nanoseconds, or -1 to wait forever.
 * @param out_frame Memory region to which the requested frame will be copied.
 * @return 0, or <0 for timeout or failure.
 */
int featbuf_consumer_wait(featbuf_t *fb, int fidx, int timeout, mfcc_t *out_frame);

/**
 * Relinquish interest in a series of frames.
 *
 * Once all consumers (defined as all objects retaining a reference to
 * @a fb) release a frame, it will no longer be available to any of
 * them.
 *
 * This can be thought of as the "up" operation for frames, where the
 * featbuf has an implicit "down" pending which removes frames from
 * circulation once they are upped enough times.
 *
 * Calling this with a non-zero timeout on the same thread which calls
 * the data processing functions below is a Bad Idea, but that should
 * be obvious.
 *
 * In order to avoid deadlock, a consumer must release all remaining
 * frames when finished with utterance processing.  This can be
 * accomplished by calling this function with @a eidx of -1.  If the
 * first frame to be released is also not known, @a sidx of 0 can also
 * be passed.
 *
 * @param fb Feature buffer.
 * @param sidx Index of first frame to be released.
 * @param eidx One past index of last frame to be released, or -1 to
 *             release all remaining frames.
 * @return 0, or <0 on error (but that is unlikely)
 */
int featbuf_consumer_release(featbuf_t *fb, int sidx, int eidx);

/**
 * Relinquish interest in a series of frames.
 *
 * This must called by a consumer thread in order to signal it is
 * finished with all utterance processing.  If the first frame to be
 * released is also not known, @a sidx of 0 can also be passed.
 *
 * @param fb Feature buffer.
 * @param sidx Index of first frame to be released.
 * @return 0, or <0 on error (but that is unlikely)
 */
int featbuf_consumer_end_utt(featbuf_t *fb, int sidx);

/**
 * Start processing for an utterance.
 *
 * @param fb Feature buffer.
 * @return 0, or <0 on error (but that is unlikely)
 */
int featbuf_producer_start_utt(featbuf_t *fb, char *uttid);

/**
 * End processing for an utterance.
 *
 * This function blocks until all consumers release the final frame of
 * the utterance.
 *
 * @param fb Feature buffer.
 * @return 0, or <0 on error (but that is unlikely)
 */
int featbuf_producer_end_utt(featbuf_t *fb);

/**
 * Shut down utterance processing.
 *
 * This function cancels any consumers waiting for an
 * utterance to start.
 *
 * @param fb Feature buffer.
 * @return 0, or <0 on error (in which case you have Problems)
 */
int featbuf_producer_shutdown(featbuf_t *fb);

/**
 * Process audio data, adding to the buffer.
 *
 * This function will always either accept all the audio passed to it,
 * or it will return an error code.
 *
 * This is not safe to call from multiple threads.
 *
 * @param fb Feature buffer.
 * @param raw Input audio data (signed 16-bit native-endian)
 * @param n_samps Number of samples in @a raw.
 * @param full_utt Does this represent an entire utterance?
 * @return 0, or <0 on error.
 */
int featbuf_producer_process_raw(featbuf_t *fb,
                                 int16 const *raw,
                                 size_t n_samps,
                                 int full_utt);

/**
 * Process acoustic feature data, adding to the buffer.
 *
 * This function will always either accept all the features passed to it,
 * or it will return an error code.
 *
 * This is not safe to call from multiple threads.
 *
 * @param fb Feature buffer.
 * @param cep Input feature data.
 * @param n_samps Number of frames in @a cep.
 * @param full_utt Does this represent an entire utterance?
 * @return 0, or <0 on error.
 */
int featbuf_producer_process_cep(featbuf_t *fb,
                                 mfcc_t **cep,
                                 size_t n_frames,
                                 int full_utt);

/**
 * Process one frame of dynamic feature data, adding to the buffer.
 *
 * This function will always either accept the features passed to it,
 * or it will return an error code.
 *
 * This is not safe to call from multiple threads.
 *
 * @param fb Feature buffer.
 * @param feat Input feature data.
 * @return 0, or <0 on error.
 */
int featbuf_producer_process_feat(featbuf_t *fb,
                                  mfcc_t **feat);

/**
 * Get the index of the first frame still being processed.
 */
int featbuf_get_window_start(featbuf_t *fb);

/**
 * Get the index of the first frame which has not yet been queued.
 */
int featbuf_get_window_end(featbuf_t *fb);

char *featbuf_uttid(featbuf_t *fb);

#endif /* __FEATBUF_H__ */
