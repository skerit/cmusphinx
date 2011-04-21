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
 * @file map_vocab.c
 * @brief Vocabulary expansion model generator.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

/* System headers. */
#include <stdio.h>

/* SphinxBase headers. */
#include <sphinxbase/pio.h>
#include <sphinxbase/err.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/filename.h>
#include <sphinxbase/byteorder.h>
#include <sphinxbase/ngram_model.h>
#include <sphinxbase/cmd_ln.h>

/* S3kr3t headerz. */
#include <multisphinx/fwdflat_search.h>
#include <multisphinx/s2_semi_mgau.h>
#include <multisphinx/cmdln_macro.h>

static const arg_t ps_args_def[] = {
    MULTISPHINX_OPTIONS,
    { "-lsn",
      ARG_STRING,
      NULL,
      "Input transcription file." },
    { "-bgdict",
      ARG_STRING,
      NULL,
      "Background dictionary file." },
    { "-prune_topn",
      ARG_INTEGER,
      NULL,
      "In absence of -bgdict, prune vocabulary by taking top-N unigrams." },
    /* Other methods will certainly folllow */
    CMDLN_EMPTY_OPTION
};

/**
 * Create synthetic senone scores using KL divergence of mixture weights.
 */
static void
kl_score_senones(acmod_t *acmod, int senid)
{
    s2_semi_mgau_t *s;
    int i, j, k;

    /* Of course this presumes semi continuous models. */
    s = (s2_semi_mgau_t *)acmod->mgau;
    /* Calculate KL divergence for all other senones. */
    memset(acmod->senone_scores, 0,
           s->n_sen * sizeof(*acmod->senone_scores));
    for (j = 0; j < s->n_feat; ++j) {
        for (k = 0; k < s->n_density; ++k) {
            float32 mixw = logmath_exp(s->lmath_8b,
                                       -s->s->mixw[j][k][senid]);
            for (i = 0; i < s->n_sen; ++i) {
                acmod->senone_scores[i]
                    += (int)(mixw * (-s->s->mixw[j][k][senid]
                                     - -s->s->mixw[j][k][i]))
                    /* fake_kl_scaled */
                    / s->n_feat;
            }
        }
    }
}

static dict_t *
pruner_topn(dict_t *fulldict, ngram_model_t *fulllm, cmd_ln_t *config, acmod_t *acmod)
{
    
    return NULL;
}

typedef struct {
    char const *arg;
    dict_t *(*func)(dict_t *fulldict, ngram_model_t *fulllm, cmd_ln_t *config, acmod_t *acmod);
} pruner_t;

static const pruner_t pruners[] = {
    { "-prune_topn", &pruner_topn },
    { NULL, NULL }
};

/**
 * Generate a pruned dictionary according to the parameters in config.
 */
static dict_t *
prune_dict(dict_t *fulldict, ngram_model_t *fulllm, cmd_ln_t *config, acmod_t *acmod)
{
    pruner_t const *p;

    for (p = pruners; p->arg; ++p) {
        if (cmd_ln_boolean_r(config, p->arg))
            return (*p->func)(fulldict, fulllm, config, acmod);
    }
    return NULL;
}

/**
 * Find exact matches between two dictionaries.
 */

/**
 * Construct a word lattice for a single unknown word.
 */

/**
 * Construct a 
 */

int
main(int argc, char *argv[])
{
#if 0
    ps_decoder_t *ps;
    cmd_ln_t *config;
    acmod_t *acmod;
    ps_search_t *fwdflat;
    dict_t *bgdict;
    char const *str;

    config = cmd_ln_parse_r(NULL, ps_args_def, argc, argv, TRUE);
    /* Get the API to initialize a bunch of stuff for us (but not the search). */
    /* FIXME: There will soon be a better way to do this. */
    cmd_ln_set_boolean_r(config, "-fwdtree", FALSE);
    cmd_ln_set_boolean_r(config, "-fwdflat", FALSE);
    cmd_ln_set_boolean_r(config, "-bestpath", FALSE);
    ps = ps_init(config);
    acmod = ps->acmod;
    /* Decide how to construct the background dictionary. */
    if ((str = cmd_ln_str_r(config, "-bgdict"))) {
        cmd_ln_t *c2 = cmd_ln_init(NULL, ps_args_def, TRUE,
                                   "-dict", str, NULL);
        bgdict = dict_init(c2, acmod->mdef);
    }
    else {
        ngram_model_t *fulllm = ngram_model_read(config,
                                                 cmd_ln_str_r(config, "-lm"),
                                                 NGRAM_AUTO, ps->lmath);
        bgdict = prune_dict(ps->dict, fulllm, config, acmod);
        ngram_model_free(fulllm);
    }

    /* Construct the fwdflat search for background. */
    fwdflat = fwdflat_search_init(config, acmod, ps->dict, ps->d2p, NULL);

    /* Compute intersection of the two vocabularies.  We first do a
     * one-to-one match on orthography (case insensitive) and then on
     * pronunciation.  Exact matches are recorded in a table.  */

    /* For non-exact matches we construct word lattices. */
    
    ps_search_free(fwdflat);
    ps_free(ps);
#endif
    return 0;
}
