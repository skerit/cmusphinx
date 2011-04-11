/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
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

#ifdef WIN32
#include <direct.h>
#endif
#include <string.h>

#include "kb.h"
#include "utt.h"
#include "srch.h"

/*
 * Begin search at bigrams of <s>, backing off to unigrams; and fillers. 
 * Update kb->lextree_next_active with the list of active lextrees.
 */
void
utt_begin(kb_t * kb)
{
    srch_utt_begin((srch_t *) kb->srch);
}

void
utt_end(kb_t * kb)
{
    srch_utt_end((srch_t *) kb->srch);
}

void
utt_decode(void *data, utt_res_t * ur, int32 sf, int32 ef, char *uttid)
{
    kb_t *kb;
    kbcore_t *kbcore;
    cmd_ln_t *config;
    int32 num_decode_frame;
    int32 total_frame;
    stat_t *st;
    srch_t *s;

    num_decode_frame = 0;
    E_INFO("Processing: %s\n", uttid);

    kb = (kb_t *) data;
    kbcore = kb->kbcore;
    config = kbcore_config(kbcore);
    kb_set_uttid(uttid, ur->uttfile, kb);
    st = kb->stat;

    /* Convert input file to cepstra if waveform input is selected */
    if (cmd_ln_boolean_r(config, "-adcin")) {
        int16 *adcdata;
        int32 nsamps = 0;

        if ((adcdata = bio_read_wavfile(cmd_ln_str_r(config, "-cepdir"),
    				        ur->uttfile,
    				        cmd_ln_str_r(config, "-cepext"),
    				        cmd_ln_int32_r(config, "-adchdr"),
    				        strcmp(cmd_ln_str_r(config, "-input_endian"), "big"),
    				        &nsamps)) == NULL) {
            E_FATAL("Cannot read file %s\n", ur->uttfile);
        }
        if (kb->mfcc) {
            ckd_free_2d((void **)kb->mfcc);
        }
        fe_start_utt(kb->fe);
        if (fe_process_utt(kb->fe, adcdata, nsamps, &kb->mfcc, &total_frame) < 0) {
            E_FATAL("MFCC calculation failed\n", ur->uttfile);
        }
        ckd_free(adcdata);
        if (total_frame > S3_MAX_FRAMES) {
            E_FATAL("Maximum number of frames (%d) exceeded\n", S3_MAX_FRAMES);
        }
        if ((total_frame = feat_s2mfc2feat_live(kbcore_fcb(kbcore),
						kb->mfcc,
						&total_frame,
						TRUE, TRUE,
						kb->feat)) < 0) {
            E_FATAL("Feature computation failed\n");
        }
    }
    else {
        /* Read mfc file and build feature vectors for entire utterance */
        if ((total_frame = feat_s2mfc2feat(kbcore_fcb(kbcore), ur->uttfile,
                                           cmd_ln_str_r(config, "-cepdir"),
                                           cmd_ln_str_r(config, "-cepext"), sf, ef,
                                           kb->feat, S3_MAX_FRAMES)) < 0) {
            E_FATAL("Cannot read file %s. Forced exit\n", ur->uttfile);
        }
    }

    /* Also need to make sure we don't set resource if it is the same. Well, this mechanism
       could be provided inside the following function. 
    */
    s = kb->srch;
    if (ur->lmname != NULL)
        srch_set_lm(s, ur->lmname);
    if (ur->regmatname != NULL)
        kb_setmllr(ur->regmatname, ur->cb2mllrname, kb);
    /* These are necessary! */
    s->uttid = kb->uttid;
    s->uttfile = kb->uttfile;

    utt_begin(kb);
    utt_decode_block(kb->feat, total_frame, &num_decode_frame, kb);
    utt_end(kb);

    st->tot_fr += st->nfr;
}


/** This function decodes a block of incoming feature vectors.
 * Feature vectors have to be computed by the calling routine.
 * The utterance level index of the last feature vector decoded
 * (before the current block) must be passed. 
 * The current status of the decode is stored in the kb structure that 
 * is passed in.
 */

void
utt_decode_block(float ***block_feat,   /* Incoming block of featurevecs */
                 int32 no_frm,  /* No. of vecs in cepblock */
                 int32 * curfrm,        /* Utterance level index of
                                           frames decoded so far */
                 kb_t * kb      /* kb structure with all model
                                   and decoder info */
    )
{

    srch_t *s;
    s = (srch_t *) kb->srch;

    /* These are necessary! */
    s->uttid = kb->uttid;
    s->uttfile = kb->uttfile;
    if (srch_utt_decode_blk(s, block_feat, no_frm, curfrm) == SRCH_FAILURE) {
        E_ERROR("srch_utt_decode_blk failed. \n");
    }
}
