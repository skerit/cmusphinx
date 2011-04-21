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
 * @file ngram_homos.c
 * @brief Merge words in an N-gram model based on a vocabulary map
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
#include <sphinxbase/garray.h>
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/bitvec.h>

/* MultiSphinx headers. */
#include <multisphinx/cmdln_macro.h>
#include <multisphinx/vocab_map.h>
#include <multisphinx/ngram_trie.h>

static const arg_t args_def[] = {
    { "-lm",
      REQARG_STRING,
      NULL,
      "Input language model file." },
    { "-vm",
      REQARG_STRING,
      NULL,
      "Input vocabulary mapping file." },
    { "-outlm",
      REQARG_STRING,
      NULL,
      "Output language model file." },
    { "-validate",
      ARG_BOOLEAN,
      "no",
      "Validate the probability distribution of the output LM." },
    CMDLN_EMPTY_OPTION
};

static int
add_weighted_successors(ngram_trie_t *lm, ngram_trie_node_t *dest,
                        ngram_trie_node_t *src, int32 weight)
{
    ngram_trie_node_t **itor;
    logmath_t *lmath;
    int32 dmarg;
    size_t nsucc, pos;

    if ((itor = ngram_trie_successors_unchecked(lm, src, &nsucc)) == NULL)
        return 0;

    lmath = ngram_trie_logmath(lm);
    E_INFOCONT("Merging (%f) [ ", logmath_exp(lmath, weight));
    ngram_trie_node_print(lm, src, err_get_logfp());

    dmarg = logmath_get_zero(lmath);
    for (pos = 0; pos < nsucc; ++pos) {
        ngram_trie_node_t *ds = itor[pos];
        int32 ds_log_prob;
        ngram_trie_node_params(lm, ds, &ds_log_prob, NULL);
        dmarg = logmath_add(lmath, dmarg, ds_log_prob);
    }
    E_INFOCONT(" ] marginal prob %g", logmath_exp(lmath, dmarg));


    for (pos = 0; pos < nsucc; ++pos) {
        ngram_trie_node_t *ss = itor[pos];
        ngram_trie_node_t *ds;
        int32 ss_log_prob;

        ngram_trie_node_params(lm, ss, &ss_log_prob, NULL);

        if ((ds = ngram_trie_successor
             (lm, dest, ngram_trie_node_word(lm, ss))) != NULL) {
            int32 ds_log_prob;
            ngram_trie_node_params(lm, ds, &ds_log_prob, NULL);
            ngram_trie_node_set_params(lm, ds,
                                       logmath_add(lmath, ds_log_prob,
                                                   weight + ss_log_prob), 0);
        }
        else {
            ngram_trie_node_set_params(lm, ss,
                                       weight + ss_log_prob, 0);
            ngram_trie_add_successor_ngram(lm, dest, ss);
        }
    }

    /* Verify marginal probability of destination. */
    dmarg = logmath_get_zero(lmath);
    for (pos = 0; pos < nsucc; ++pos) {
        ngram_trie_node_t *ds = itor[pos];
        int32 ds_log_prob;

        ngram_trie_node_params(lm, ds, &ds_log_prob, NULL);
        dmarg = logmath_add(lmath, dmarg, ds_log_prob);
    }
    E_INFOCONT(" dest marginal prob %d = %g\n",
               dmarg, logmath_exp(lmath, dmarg));

    return 0;
}

typedef struct fast_vmap_s {
    int32 pseudo_wid;
    int32 n_wids;
    bitvec_t *wids;
} fast_vmap_t;

static int32
find_word_succ(ngram_trie_t *lm, ngram_trie_node_t *node,
               garray_t *word_succ, fast_vmap_t *vmap)
{
    int32 succ_prob;
    dict_t *dict;
    logmath_t *lmath;
    ngram_trie_node_t **ni;
    size_t nsucc, pos;

    dict = ngram_trie_dict(lm);
    lmath = ngram_trie_logmath(lm);
    succ_prob = logmath_get_zero(lmath);
    if ((ni = ngram_trie_successors_unchecked(lm, node, &nsucc)) == NULL)
        return succ_prob;

    for (pos = 0; pos < nsucc; ++pos) {
        ngram_trie_node_t *succ = ni[pos];
        int32 basewid = ngram_trie_node_word(lm, succ);
        if (bitvec_is_set(vmap->wids, basewid)) {
            int32 log_prob, log_bowt;
            i32p_t sent;
            ngram_trie_node_params(lm, succ, &log_prob, &log_bowt);
            sent.a = basewid;
            sent.b = log_prob;
            garray_append(word_succ, &sent);
            succ_prob = logmath_add(lmath, succ_prob, log_prob);
        }
    }
    return succ_prob;
}

static int
merge_successors(ngram_trie_t *lm, ngram_trie_node_t *node,
                 garray_t *word_succ, int32 succ_prob, int32 pseudo_wid)
{
    ngram_trie_node_t *pseudo_succ;
    int i;
    pseudo_succ = ngram_trie_node_alloc(lm);
    ngram_trie_node_set_word(lm, pseudo_succ, pseudo_wid);
    ngram_trie_node_set_params(lm, pseudo_succ, succ_prob, 0);
    for (i = 0; i < garray_size(word_succ); ++i) {
        i32p_t sent = garray_ent(word_succ, i32p_t, i);
        ngram_trie_node_t *succ;
        int32 weight;
        succ = ngram_trie_successor(lm, node, sent.a);
        weight = sent.b - succ_prob;
        if (add_weighted_successors(lm, pseudo_succ,
                                    succ, weight) < 0) {
            garray_free(word_succ);
            return -1;
        }
        if (ngram_trie_delete_successor(lm, node, sent.a) < 0) {
            E_ERROR("Failed to delete successor\n");
            return -1;
        }
    }
    return ngram_trie_add_successor_ngram(lm, node, pseudo_succ);
}


static int
merge_homos(ngram_trie_t *lm, ngram_trie_node_t *node, garray_t *vmaps)
{
    ngram_trie_iter_t *ni;
    logmath_t *lmath;
    garray_t *word_succ;
    dict_t *dict;
    size_t pos, nmaps;
    fast_vmap_t *vmaps_arr;

    /* Stop condition: This is a leaf node, nothing to do. */
    if ((ni = ngram_trie_successors(lm, node)) == NULL)
        return 0;

    /* Traverse the N-Gram trie in pre-order. */
    for (; ni; ni = ngram_trie_iter_next(ni))
        if (merge_homos(lm, ngram_trie_iter_get(ni), vmaps) < 0)
            return -1;

    /* Merge within-distribution probabilities. */
    dict = ngram_trie_dict(lm);
    lmath = ngram_trie_logmath(lm);
    word_succ = garray_init(0, sizeof(i32p_t));
    garray_set_cmp(word_succ, garray_cmp_i32p_first, NULL);
    nmaps = garray_size(vmaps);
    vmaps_arr = garray_ptr(vmaps, fast_vmap_t, 0);
    for (pos = 0; pos < nmaps; ++pos) {
        fast_vmap_t *vmap = vmaps_arr + pos;
        int32 succ_prob;

        /* Collect base words to be merged (if any). */
        garray_reset(word_succ);
        if (vmap->n_wids == 0)
            continue;
        /* Calculate normalizer for weights along the way. */
        succ_prob = find_word_succ(lm, node, word_succ, vmap);
        if (garray_size(word_succ) == 0)
            continue;
        else if (garray_size(word_succ) == 1) {
            i32p_t sent = garray_ent(word_succ, i32p_t, 0);
            ngram_trie_node_t *pseudo_succ;
            pseudo_succ = ngram_trie_successor(lm, node, sent.a);
            ngram_trie_rename_successor(lm, node, pseudo_succ, vmap->pseudo_wid);
        }
        else {
            if (merge_successors(lm, node, word_succ, succ_prob, vmap->pseudo_wid) < 0)
                return -1;
        }
    }
    garray_free(word_succ);
    return 0;
}

static int
recalc_bowts(ngram_trie_t *lm)
{
    ngram_trie_iter_t *ng;
    int n;

    for (n = 1; n < ngram_trie_n(lm); ++n) {
        for (ng = ngram_trie_ngrams(lm, n); ng;
             ng = ngram_trie_iter_next(ng)) {
            int32 old_bowt, log_bowt, log_prob;
            ngram_trie_node_t *node = ngram_trie_iter_get(ng);

            ngram_trie_node_params(lm, node, &log_prob, &old_bowt);
	    /* ngram_trie_node_print(lm, node, stdout); */
            log_bowt = ngram_trie_calc_bowt(lm, node);
	    /* printf(": %d %d\n", old_bowt, log_bowt); */
            ngram_trie_node_set_params(lm, node, log_prob, log_bowt);
        }
    }
    return 0;
}

static int
validate(ngram_trie_t *lm)
{
    ngram_trie_iter_t *ng;
    int n;

    for (n = 1; n < ngram_trie_n(lm); ++n) {
        for (ng = ngram_trie_ngrams(lm, n); ng;
             ng = ngram_trie_iter_next(ng)) {
            ngram_trie_node_t *node = ngram_trie_iter_get(ng);
            if (!ngram_trie_node_validate(lm, node))
                return FALSE;
        }
    }
    return TRUE;
}

int
main(int argc, char *argv[])
{
    cmd_ln_t *config;
    vocab_map_t *vm;
    vocab_map_iter_t *vi;
    garray_t *vmaps;
    ngram_trie_t *lm;
    logmath_t *lmath;
    dict_t *dict, *vdict;
    FILE *fh;
    int32 ispipe;

    config = cmd_ln_parse_r(NULL, args_def, argc, argv, TRUE);
    if (config == NULL)
        return 1;

    lmath = logmath_init(1.0003, 0, FALSE);
    lm = ngram_trie_init(NULL, lmath);
    if ((fh = fopen_comp(cmd_ln_str_r(config, "-lm"), "r", &ispipe)) == NULL)
        E_FATAL_SYSTEM("Failed to open %s", cmd_ln_str_r(config, "-lm"));
    if (ngram_trie_read_arpa(lm, fh) < 0)
        return 1;
    fclose(fh);

    vm = vocab_map_init(NULL);
    if ((fh = fopen(cmd_ln_str_r(config, "-vm"), "r")) == NULL)
        E_FATAL_SYSTEM("Failed to open %s", cmd_ln_str_r(config, "-vm"));
    if (vocab_map_read(vm, fh) < 0)
        return 1;
    fclose(fh);

    /* Add new pseudo-words to the LM dictionary. */
    dict = ngram_trie_dict(lm);
    vdict = vocab_map_dict(vm);
    for (vi = vocab_map_mappings(vm); vi;
         vi = vocab_map_iter_next(vi)) {
        char const *pseudo_word;
        fast_vmap_t vmap;
        int32 const *wids;

        wids = vocab_map_iter_get(vi, &vmap.pseudo_wid, &vmap.n_wids);
        pseudo_word = dict_wordstr(vdict, vmap.pseudo_wid);
        /* FIXME: No pronunciation for now. */
        dict_add_word(dict, pseudo_word, NULL, 0);
    }
    /* Now re-read the vocab map with this updated dictionary. */
    vocab_map_free(vm);
    vm = vocab_map_init(dict);
    if ((fh = fopen(cmd_ln_str_r(config, "-vm"), "r")) == NULL)
        E_FATAL_SYSTEM("Failed to open %s", cmd_ln_str_r(config, "-vm"));
    if (vocab_map_read(vm, fh) < 0)
        return 1;
    fclose(fh);

    /* Construct fast word mappings. */
    vmaps = garray_init(0, sizeof(fast_vmap_t));
    for (vi = vocab_map_mappings(vm); vi;
         vi = vocab_map_iter_next(vi)) {
        fast_vmap_t vmap;
        int32 const *wids;
        int32 i;

        wids = vocab_map_iter_get(vi, &vmap.pseudo_wid, &vmap.n_wids);
        vmap.wids = bitvec_alloc(dict_size(dict));
        for (i = 0; i < vmap.n_wids; ++i)
            bitvec_set(vmap.wids, wids[i]);
        garray_append(vmaps, &vmap);
    }
    if (merge_homos(lm, ngram_trie_root(lm), vmaps) < 0)
        return 1;
    if (recalc_bowts(lm) < 0)
        return 1;

    if ((fh = fopen_comp(cmd_ln_str_r(config, "-outlm"), "w", &ispipe)) == NULL)
        E_FATAL_SYSTEM("Failed to open %s", cmd_ln_str_r(config, "-outlm"));
    if (ngram_trie_write_arpa(lm, fh) < 0)
        return 1;
    if (fclose(fh) < 0)
        E_FATAL_SYSTEM("Failed to complete writing ARPA file");

    if (cmd_ln_boolean_r(config, "-validate"))
        if (!validate(lm))
            return 1;

    vocab_map_free(vm);
    ngram_trie_free(lm);
    cmd_ln_free_r(config);
    logmath_free(lmath);

    return 0;
}
