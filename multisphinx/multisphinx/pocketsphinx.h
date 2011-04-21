/* -*- c-basic-offset:4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2010 Carnegie Mellon University.  All rights
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
 * @file pocketsphinx.h PocketSphinx API compatibility
 */

#ifndef __POCKETSPHINX_H__
#define __POCKETSPHINX_H__

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/* System headers we need. */
#include <stdio.h>

/* SphinxBase headers we need. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/fe.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/ngram_model.h>

/* PocketSphinx headers (not many of them!) */
#include <multisphinx/win32_export.h>
#include <multisphinx/cmdln_macro.h>

/**
 * PocketSphinx speech recognizer object.
 */
typedef struct ps_decoder_s ps_decoder_t;

/**
 * PocketSphinx segmentation iterator object.
 */
typedef struct ps_seg_s ps_seg_t;

/**
 * Initialize the decoder from a configuration object.
 *
 * @note The decoder retains ownership of the pointer
 * <code>config</code>, so you must not attempt to free it manually.
 * If you wish to reuse it elsewhere, call cmd_ln_retain() on it.
 *
 * @param config a command-line structure, as created by
 * cmd_ln_parse_r() or cmd_ln_parse_file_r().
 */
POCKETSPHINX_EXPORT
ps_decoder_t *ps_init(cmd_ln_t *config);

/**
 * Returns the argument definitions used in ps_init().
 *
 * This is here to avoid exporting global data, which is problematic
 * on Win32 and Symbian (and possibly other platforms).
 */
POCKETSPHINX_EXPORT
arg_t const *ps_args(void);

/**
 * Set up default arguments in a command-line object.
 *
 * This function initializes a command-line with default acoustic and
 * language models if not already set.  It also finds the relevant
 * acoustic model files if only -hmm is provided.
 */
POCKETSPHINX_EXPORT
void ps_init_defaults(cmd_ln_t *config);

/**
 * Retain a pointer to the decoder.
 *
 * This increments the reference count on the decoder, allowing it to
 * be shared between multiple parent objects.  In general you will not
 * need to use this function, ever.  It is mainly here for the
 * convenience of scripting language bindings.
 *
 * @return pointer to retained decoder.
 */
POCKETSPHINX_EXPORT
ps_decoder_t *ps_retain(ps_decoder_t *ps);

/**
 * Finalize the decoder.
 *
 * This releases all resources associated with the decoder, including
 * any language models or grammars which have been added to it, and
 * the initial configuration object passed to ps_init().
 *
 * @param ps Decoder to be freed.
 * @return New reference count (0 if freed).
 */
POCKETSPHINX_EXPORT
int ps_free(ps_decoder_t *ps);

/**
 * Get the configuration object for this decoder.
 *
 * @return The configuration object for this decoder.  The decoder
 *         retains ownership of this pointer, so you should not
 *         attempt to free it manually.  Use cmd_ln_retain() if you
 *         wish to reuse it elsewhere.
 */
POCKETSPHINX_EXPORT
cmd_ln_t *ps_get_config(ps_decoder_t *ps);

/**
 * Get the log-math computation object for this decoder.
 *
 * @return The log-math object for this decoder.  The decoder retains
 *         ownership of this pointer, so you should not attempt to
 *         free it manually.  Use logmath_retain() if you wish to
 *         reuse it elsewhere.
 */
POCKETSPHINX_EXPORT
logmath_t *ps_get_logmath(ps_decoder_t *ps);

/**
 * Get the feature extraction object for this decoder.
 *
 * @return The feature extraction object for this decoder.  The
 *         decoder retains ownership of this pointer, so you should
 *         not attempt to free it manually.  Use fe_retain() if you
 *         wish to reuse it elsewhere.
 */
POCKETSPHINX_EXPORT
fe_t *ps_get_fe(ps_decoder_t *ps);

/**
 * Get the dynamic feature computation object for this decoder.
 *
 * @return The dynamic feature computation object for this decoder.  The
 *         decoder retains ownership of this pointer, so you should
 *         not attempt to free it manually.  Use feat_retain() if you
 *         wish to reuse it elsewhere.
 */
POCKETSPHINX_EXPORT
feat_t *ps_get_feat(ps_decoder_t *ps);

/**
 * Decode a raw audio stream.
 *
 * No headers are recognized in this files.  The configuration
 * parameters <tt>-samprate</tt> and <tt>-input_endian</tt> are used
 * to determine the sampling rate and endianness of the stream,
 * respectively.  Audio is always assumed to be 16-bit signed PCM.
 *
 * @param ps Decoder.
 * @param rawfh Previously opened file stream.
 * @param uttid Utterance ID (or NULL to generate automatically).
 * @param maxsamps Maximum number of samples to read from rawfh, or -1
 *                 to read until end-of-file.
 * @return Number of samples of audio.
 */
POCKETSPHINX_EXPORT
int ps_decode_raw(ps_decoder_t *ps, FILE *rawfh,
                  char const *uttid, long maxsamps);

/**
 * Start utterance processing.
 *
 * This function should be called before any utterance data is passed
 * to the decoder.  It marks the start of a new utterance and
 * reinitializes internal data structures.
 *
 * @param ps Decoder to be started.
 * @param uttid String uniquely identifying this utterance.  If NULL,
 *              one will be created.
 * @return 0 for success, <0 on error.
 */
POCKETSPHINX_EXPORT
int ps_start_utt(ps_decoder_t *ps, char const *uttid);

/**
 * Get current utterance ID.
 *
 * @param ps Decoder to query.
 * @return Read-only string of the current utterance ID.  This is
 * valid only until the beginning of the next utterance.
 */
POCKETSPHINX_EXPORT
char const *ps_get_uttid(ps_decoder_t *ps);
/**
 * Decode raw audio data.
 *
 * @param ps Decoder.
 * @param no_search This option does nothing.  This function never blocks.
 * @param full_utt If non-zero, this block of data is a full utterance
 *                 worth of data.  This may allow the recognizer to
 *                 produce more accurate results.
 * @return 0, or <0 for error.
 */
POCKETSPHINX_EXPORT
int ps_process_raw(ps_decoder_t *ps,
                   int16 const *data,
                   size_t n_samples,
                   int no_search,
                   int full_utt);

/**
 * Decode acoustic feature data.
 *
 * @param ps Decoder.
 * @param no_search This option does nothing.  This function never blocks.
 * @param full_utt If non-zero, this block of data is a full utterance
 *                 worth of data.  This may allow the recognizer to
 *                 produce more accurate results.
 * @return 0, or <0 for error.
 */
POCKETSPHINX_EXPORT
int ps_process_cep(ps_decoder_t *ps,
                   mfcc_t **data,
                   int n_frames,
                   int no_search,
                   int full_utt);

/**
 * Get the number of frames of data searched.
 *
 * Note that there is a delay between this and the number of frames of
 * audio which have been input to the system.  This is due to the fact
 * that acoustic features are computed using a sliding window of
 * audio, and dynamic features are computed over a sliding window of
 * acoustic features.
 *
 * @param ps Decoder.
 * @return Number of frames of speech data which have been recognized
 * so far.
 */
POCKETSPHINX_EXPORT
int ps_get_n_frames(ps_decoder_t *ps);

/**
 * End utterance processing.
 *
 * @param ps Decoder.
 * @return 0 for success, <0 on error
 */
POCKETSPHINX_EXPORT
int ps_end_utt(ps_decoder_t *ps);

/**
 * Get hypothesis string and path score.
 *
 * @param ps Decoder.
 * @param out_best_score Output: path score corresponding to returned string.
 * @param out_uttid Output: utterance ID for this utterance.
 * @return String containing best hypothesis at this point in
 *         decoding.  NULL if no hypothesis is available.
 */
POCKETSPHINX_EXPORT
char const *ps_get_hyp(ps_decoder_t *ps, int32 *out_best_score,
                       char const **out_uttid);

/**
 * Get posterior probability.
 *
 * @note Unless the -bestpath option is enabled, this function will
 * always return zero (corresponding to a posterior probability of
 * 1.0).  Even if -bestpath is enabled, it will also return zero when
 * called on a partial result.  Ongoing research into effective
 * confidence annotation for partial hypotheses may result in these
 * restrictions being lifted in future versions.
 *
 * @param ps Decoder.
 * @param out_uttid Output: utterance ID for this utterance.
 * @return Posterior probability of the best hypothesis.
 */
POCKETSPHINX_EXPORT
int32 ps_get_prob(ps_decoder_t *ps, char const **out_uttid);

/**
 * Get an iterator over the word segmentation for the best hypothesis.
 *
 * @param ps Decoder.
 * @param out_best_score Output: path score corresponding to hypothesis.
 * @return Iterator over the best hypothesis at this point in
 *         decoding.  NULL if no hypothesis is available.
 */
POCKETSPHINX_EXPORT
ps_seg_t *ps_seg_iter(ps_decoder_t *ps, int32 *out_best_score);

/**
 * Get the next segment in a word segmentation.
 *
 * @param seg Segment iterator.
 * @return Updated iterator with the next segment.  NULL at end of
 *         utterance (the iterator will be freed in this case).
 */
POCKETSPHINX_EXPORT
ps_seg_t *ps_seg_next(ps_seg_t *seg);

/**
 * Get word string from a segmentation iterator.
 *
 * @param seg Segment iterator.
 * @return Read-only string giving string name of this segment.  This
 * is only valid until the next call to ps_seg_next().
 */
POCKETSPHINX_EXPORT
char const *ps_seg_word(ps_seg_t *seg);

/**
 * Get inclusive start and end frames from a segmentation iterator.
 *
 * @note These frame numbers are inclusive, i.e. the end frame refers
 * to the last frame in which the given word or other segment was
 * active.  Therefore, the actual duration is *out_ef - *out_sf + 1.
 *
 * @param seg Segment iterator.
 * @param out_sf Output: First frame index in segment.
 * @param out_sf Output: Last frame index in segment.
 */
POCKETSPHINX_EXPORT
void ps_seg_frames(ps_seg_t *seg, int *out_sf, int *out_ef);

/**
 * Get language, acoustic, and posterior probabilities from a
 * segmentation iterator.
 *
 * @note Unless the -bestpath option is enabled, this function will
 * always return zero (corresponding to a posterior probability of
 * 1.0).  Even if -bestpath is enabled, it will also return zero when
 * called on a partial result.  Ongoing research into effective
 * confidence annotation for partial hypotheses may result in these
 * restrictions being lifted in future versions.
 *
 * @param out_ascr Output: acoustic model score for this segment.
 * @param out_lscr Output: language model score for this segment.
 * @param out_lback Output: language model backoff mode for this
 *                  segment (i.e. the number of words used in
 *                  calculating lscr).  This field is, of course, only
 *                  meaningful for N-Gram models.
 * @return Log posterior probability of current segment.  Log is
 *         expressed in the log-base used in the decoder.  To convert
 *         to linear floating-point, use logmath_exp(ps_get_logmath(),
 *         pprob).
 */
POCKETSPHINX_EXPORT
int32 ps_seg_prob(ps_seg_t *seg, int32 *out_ascr, int32 *out_lscr, int32 *out_lback);

/**
 * Finish iterating over a word segmentation early, freeing resources.
 */
POCKETSPHINX_EXPORT
void ps_seg_free(ps_seg_t *seg);

/**
 * Get performance information for the current utterance.
 *
 * @param ps Decoder.
 * @param out_nspeech Output: Number of seconds of speech.
 * @param out_ncpu    Output: Number of seconds of CPU time used.
 * @param out_nwall   Output: Number of seconds of wall time used.
 */
POCKETSPHINX_EXPORT
void ps_get_utt_time(ps_decoder_t *ps, double *out_nspeech,
                     double *out_ncpu, double *out_nwall);

/**
 * Get overall performance information.
 *
 * @param ps Decoder.
 * @param out_nspeech Output: Number of seconds of speech.
 * @param out_ncpu    Output: Number of seconds of CPU time used.
 * @param out_nwall   Output: Number of seconds of wall time used.
 */
POCKETSPHINX_EXPORT
void ps_get_all_time(ps_decoder_t *ps, double *out_nspeech,
                     double *out_ncpu, double *out_nwall);

/**
 * @mainpage PocketSphinx API Documentation
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 * @version 0.8
 * @date October, 2010
 *
 * @section intro_sec Introduction
 *
 * This is the API documentation for the PocketSphinx speech
 * recognition engine.  The main API calls are documented in
 * <pocketsphinx.h>.
 *
 * Please note that the API has changed since the previous revision,
 * such that it may no longer be binary-compatible with your program.
 * The behaviour of some API calls is different as well.
 *
 * The main difference is that all calls which involve processing data
 * are now non-blocking.  Therefore, the @a no_search parameter no
 * longer has any effect.  In addition, this means that ps_get_hyp()
 * and ps_seg_iter() return the most recent result available, which
 * may not correspond to the last block of audio or features
 * submitted, depending on the speed of your platform.
 *
 * As previously, ps_end_utt() is blocking (although this wasn't
 * explicitly stated before) and will block until the final hypothesis
 * is available.  Due to the "anytime" nature of PocketSphinx
 * recognition, a full hypothesis may be available much earlier.
 * Therefore, a new function, ps_wait(), has been added, which blocks
 * until a certain timepoint in the input has a result available, with
 * an option to wait until the final result for that timepoint is
 * available.  Additionally, this function can be used to wait for all
 * data submitted so far to be processed.
 *
 * The internal handling of dictionaries and language models has
 * changed to allow dynamic vocabularies and domains to be handled in
 * a more accurate and principled fashion.  The ps_vocab_t interface
 * is provided to allow vocabulary manipulation at runtime. (FIXME: at
 * least, it will be...)
 */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __POCKETSPHINX_H__ */
