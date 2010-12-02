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

#ifndef __VOCAB_MAPPER_H__
#define __VOCAB_MAPPER_H__

#include "vocab_map.h"

typedef struct vocab_mapper_s vocab_mapper_t;

/**
 * Initialize a vocabulary mapper.
 *
 * Creates a new vocabulary mapper from the specified target language
 * model and dictionary.
 *
 * @param config Configuration parameters.  Recognized parameters are:
 *  - @a hmm Acoustic model to use for mapping
 *  - @a prune_topn Prune background model to this number of word types

 * @param targdict Target dictionary, i.e. the dictionary which will
 *                 ultimately be recognized.
 * @param targlm Target language model, used in final recognition.
 * @return Newly created vocabulary mapper.
 */
vocab_mapper_t *vocab_mapper_init(cmd_ln_t *config,
				  dict_t *targdict,
				  ngram_model_t *targlm);

/**
 * Retain a pointer to a vocabulary mapper.
 */
vocab_mapper_t *vocab_mapper_retain(vocab_mapper_t *vm);

/**
 * Release a pointer to a vocabulary mapper.
 */
int vocab_mapper_free(vocab_mapper_t *vm);

/**
 * Set explicit background language model.
 */
int vocab_mapper_set_bglm(vocab_mapper_t *vm, ngram_model_t *bglm);

/**
 * Set explicit background dictionary.
 */
int vocab_mapper_set_bgdict(vocab_mapper_t *vm, dict_t *bgdict);

/**
 * Generate background language model and dictionary by pruning based
 * on unigram probability.
 */
int vocab_mapper_prune_unigram(vocab_mapper_t *vm);

/**
 * Get background language model.
 */
ngram_model_t *vocab_mapper_get_bglm(vocab_mapper_t *vm);

/**
 * Get background dictionary.
 */
dict_t *vocab_mapper_get_bgdict(vocab_mapper_t *vm);

/**
 * Get vocabulary map.
 */
vocab_map_t *vocab_mapper_get_vocab_map(vocab_mapper_t *vm);

#endif /* __VOCAB_MAPPER_H__ */
