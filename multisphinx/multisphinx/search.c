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

#include "search_internal.h"

void
search_base_init(search_t *search, searchfuncs_t *vt,
            cmd_ln_t *config, acmod_t *acmod, dict2pid_t *d2p)
{
    ptmr_init(&search->t);
    search->vt = vt;
    search->config = cmd_ln_retain(config);
    if (acmod)
        search->acmod = acmod_retain(acmod);
    if (d2p) {
        search->d2p = dict2pid_retain(d2p);
        search->dict = dict_retain(d2p->dict);
        search->start_wid = dict_startwid(d2p->dict);
        search->finish_wid = dict_finishwid(d2p->dict);
        search->silence_wid = dict_silwid(d2p->dict);
        search->n_words = dict_size(d2p->dict);
    }
    else {
        search->d2p = NULL;
        search->dict = NULL;
        search->start_wid = search->finish_wid = search->silence_wid = -1;
        search->n_words = 0;
    }
    search->mtx = sbmtx_init();
}

int
search_free(search_t *search)
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

char const *
search_name(search_t *search)
{
    return search->vt->name;
}


arc_buffer_t *
search_link(search_t *from, search_t *to,
               char const *name, int keep_scores)
{
    arc_buffer_t *ab;

    if (search_bptbl(from) == NULL)
        return NULL;
    ab = arc_buffer_init(name, search_bptbl(from),
                         search_lmset(from), keep_scores);
    search_output_arcs(from) = ab;
    search_input_arcs(to) = arc_buffer_retain(ab);

    return ab;
}

static int
search_main(sbthread_t *thr)
{
    search_t *search = sbthread_arg(thr);
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
search_run(search_t *search)
{
    search->thr = sbthread_start(NULL, search_main, search);
    return search->thr;
}

int
search_call_event(search_t *search, int event, int frame)
{
    search_event_t evt;
    if (search->cb != NULL) {
        evt.event = event;
        evt.frame = frame;
        return (*search->cb)(search, &evt, search->cb_data);
    }
    else {
        return 0;
    }
}

void
search_set_cb(search_t *search, search_cb_func cb, void *udata)
{
    search->cb = cb;
    search->cb_data = udata;
}


int
search_wait(search_t *search)
{
    return sbthread_wait(search->thr);
}

char const *
search_hyp(search_t *search, int32 *out_score)
{
    /* Search module has to be responsible for locking here... */
    return (*search->vt->hyp)(search, out_score);
}

seg_iter_t *
search_seg_iter(search_t *search, int32 *out_score)
{
    return (*search->vt->seg_iter)(search, out_score);
}

seg_iter_t *
seg_iter_next(seg_iter_t *itor)
{
    return (*itor->vt->seg_next)(itor);
}

void
seg_iter_free(seg_iter_t *itor)
{
    (*itor->vt->seg_free)(itor);
}

char const *
seg_iter_word(seg_iter_t *itor)
{
    return dict_wordstr(search_dict(itor->search), itor->wid);
}

int32
seg_iter_wid(seg_iter_t *itor)
{
    return itor->wid;
}

void
seg_iter_times(seg_iter_t *itor, int *out_sf, int *out_ef)
{
    if (out_sf) *out_sf = itor->sf;
    if (out_ef) *out_ef = itor->ef;
}

bptbl_t *
search_bptbl(search_t *search)
{
    if (search->vt->bptbl == NULL)
        return NULL;
    return (*search->vt->bptbl)(search);
}

ngram_model_t *
search_lmset(search_t *search)
{
    if (search->vt->lmset == NULL)
        return NULL;
    return (*search->vt->lmset)(search);
}

cmd_ln_t *
search_config(search_t *search)
{
    return search->config;
}

char const *
search_uttid(search_t *search)
{
    return search->uttid;
}
