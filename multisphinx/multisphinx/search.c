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
 * @file ps_search.c Search algorithm class
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include <multisphinx/search.h>
#include <multisphinx/arc_buffer.h>

void
ps_search_init(ps_search_t *search, ps_searchfuncs_t *vt,
               cmd_ln_t *config, acmod_t *acmod, dict_t *dict,
               dict2pid_t *d2p)
{
    ptmr_init(&search->t);
    search->vt = vt;
    search->config = cmd_ln_retain(config);
    if (acmod)
        search->acmod = acmod_retain(acmod);
    if (d2p)
        search->d2p = dict2pid_retain(d2p);
    else
        search->d2p = NULL;
    if (dict) {
        search->dict = dict_retain(dict);
        search->start_wid = dict_startwid(dict);
        search->finish_wid = dict_finishwid(dict);
        search->silence_wid = dict_silwid(dict);
        search->n_words = dict_size(dict);
    }
    else {
        search->dict = NULL;
        search->start_wid = search->finish_wid = search->silence_wid = -1;
        search->n_words = 0;
    }
    search->mtx = sbmtx_init();
}

int
ps_search_free(ps_search_t *search)
{
    if (search == NULL)
        return 0;

    ptmr_stop(&search->t);
    E_INFO("TOTAL %s %.3f wall %.3f xRT\n",
           search->vt->name, search->t.t_tot_elapsed,
           search->t.t_tot_elapsed / search->total_frames
           * cmd_ln_int32_r(search->config, "-frate"));
    /* Call the search free function. */
    (*search->vt->free)(search);
    /* Clean up common stuff. */
    arc_buffer_free(search->input_arcs);
    arc_buffer_free(search->output_arcs);
    cmd_ln_free_r(search->config);
    acmod_free(search->acmod);
    dict_free(search->dict);
    dict2pid_free(search->d2p);
    ckd_free(search->hyp_str);
    sbthread_free(search->thr);
    sbmtx_free(search->mtx);
    ckd_free(search);
    return 0;
}

arc_buffer_t *
ps_search_link(ps_search_t *from, ps_search_t *to,
               char const *name, int keep_scores)
{
    arc_buffer_t *ab;

    if (ps_search_bptbl(from) == NULL)
        return NULL;
    ab = arc_buffer_init(name, ps_search_bptbl(from),
                         ps_search_lmset(from), keep_scores);
    ps_search_output_arcs(from) = ab;
    ps_search_input_arcs(to) = arc_buffer_retain(ab);

    return ab;
}

static int
ps_search_main(sbthread_t *thr)
{
    ps_search_t *search = sbthread_arg(thr);
    int rv = 0;

    while (rv >= 0) {
        ptmr_reset(&search->t);
        if ((rv = (*search->vt->decode)(search)) < 0) {
            E_INFO("%s canceled\n", search->vt->name);
        }
    }
    return 0;
}

sbthread_t *
ps_search_run(ps_search_t *search)
{
    search->thr = sbthread_start(NULL, ps_search_main, search);
    return search->thr;
}

int
ps_search_wait(ps_search_t *search)
{
    return sbthread_wait(search->thr);
}

char const *
ps_search_hyp(ps_search_t *search, int32 *out_score)
{
    /* Search module has to be responsible for locking here... */
    return (*search->vt->hyp)(search, out_score);
}

char const *
ps_search_splice(ps_search_t **searches, int nsearches,
                 int32 *out_score)
{
    /* How to do this...  Requires some cooperation from the searches.
     * Basically instead of ps_search_hyp() they should always return
     * segmentations.  Then we pass these off to some other function
     * in here which concatenates the segmentations, and we resolve
     * overlaps according to some heuristic to be determined.
     */
    return NULL;
}

ps_seg_t *
ps_search_seg_iter(ps_search_t *search, int32 *out_score)
{
    return (*search->vt->seg_iter)(search, out_score);
}

bptbl_t *
ps_search_bptbl(ps_search_t *search)
{
    if (search->vt->bptbl == NULL)
        return NULL;
    return (*search->vt->bptbl)(search);
}

ngram_model_t *
ps_search_lmset(ps_search_t *search)
{
    if (search->vt->lmset == NULL)
        return NULL;
    return (*search->vt->lmset)(search);
}
