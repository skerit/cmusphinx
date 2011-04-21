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
 * @file vocab_mapper.h
 * @brief Vocabulary mapping
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include "vocab_mapper.h"

struct vocab_mapper_s {
    int refcount;
    dict_t *bgdict;         /**< Background dictionary */
    ngram_model_t *bglm;    /**< Background language model */
    dict_t *targdict;       /**< Target dictionary */
    ngram_model_t *targlm;  /**< Target language model */
    search_t *fwdflat;   /**< Search used to generate mappings */
    acmod_t *acmod;         /**< Acoustic model used in search */
} vocab_mapper_t;

vocab_mapper_t *
vocab_mapper_init(cmd_ln_t *config,
                  dict_t *targdict,
                  ngram_model_t *targlm)
{
}

vocab_mapper_t *
vocab_mapper_retain(vocab_mapper_t *vm)
{
}

int
vocab_mapper_free(vocab_mapper_t *vm)
{
}

int
vocab_mapper_set_bglm(vocab_mapper_t *vm, ngram_model_t *bglm)
{
}

int
vocab_mapper_set_bgdict(vocab_mapper_t *vm, dict_t *bgdict)
{
}

int
vocab_mapper_prune_unigram(vocab_mapper_t *vm)
{
}

ngram_model_t *
vocab_mapper_get_bglm(vocab_mapper_t *vm)
{
}

vocab_map_t *
vocab_mapper_get_vocab_map(vocab_mapper_t *vm)
{
}
