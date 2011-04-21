/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2001 Carnegie Mellon University.  All rights
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

/* System headers. */
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

/* SphinxBase headers. */
#include <sphinxbase/pio.h>
#include <sphinxbase/err.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/filename.h>
#include <sphinxbase/byteorder.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/feat.h>

/* MultiSphinx headers. */
#include <multisphinx/cmdln_macro.h>
#include <multisphinx/dict.h>
#include <multisphinx/search.h>
#include <multisphinx/search_factory.h>
#include <multisphinx/alignment.h>
#include <multisphinx/vocab_map.h>
#include <multisphinx/fwdflat_search.h>

typedef struct batch_decoder_s {
    search_factory_t *sf;
    cmd_ln_t *config;
    search_t *fwdtree;
    search_t *fwdflat;
    search_t *latgen;

    struct timeval utt_start;

    FILE *ctlfh;
    FILE *alignfh;
    FILE *hypfh;

    hash_table_t *hypfiles;
} batch_decoder_t;

/**
 * Command-line argument definitions for batch processing.
 */
static const arg_t ms_args_def[] = { MULTISPHINX_OPTIONS,
/* Control file. */
{ "-ctl",
ARG_STRING,
NULL,
"Control file listing utterances to be processed"},
{"-ctloffset",
ARG_INT32,
"0",
"No. of utterances at the beginning of -ctl file to be skipped"},
{"-ctlcount",
ARG_INT32,
"-1",
"No. of utterances to be processed (after skipping -ctloffset entries)"},
{"-ctlincr",
ARG_INT32,
"1",
"Do every Nth line in the control file"},
/* Alignment file with times for input words */
{"-align",
ARG_STRING,
NULL,
"Alignment file with input word times for latency measurement"},
/* Prefix for hypothesis files */
{"-hypprefix",
ARG_STRING,
NULL,
"Prefix for partial hypothesis files (will be suffixed with .{pass}.hyp"},
{"-hyp",
ARG_STRING,
NULL,
"Final hypothesis file."},

/* Input file types and locations. */
{"-adcin",
ARG_BOOLEAN,
"no",
"Input is raw audio data"},
{"-adchdr",
ARG_INT32,
"0",
"Size of audio file header in bytes (headers are ignored)"},
{"-cepdir",
ARG_STRING,
NULL,
"Input files directory (prefixed to filespecs in control file)"},
{"-cepext",
ARG_STRING,
".mfc",
"Input files extension (suffixed to filespecs in control file)"},

/* Terminator */
CMDLN_EMPTY_OPTION
};

static double get_time_delta(batch_decoder_t *bd)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return ((double) tv.tv_sec + (double) tv.tv_usec / 1000000)
            - ((double) bd->utt_start.tv_sec + (double) bd->utt_start.tv_usec
                    / 1000000);
}

static int batch_decoder_decode_adc(batch_decoder_t *bd, FILE *infh, int sf,
        int ef, alignment_t *al)
{
    featbuf_t *fb = search_factory_featbuf(bd->sf);
    float32 samprate = cmd_ln_float32_r(bd->config, "-samprate");
    int32 frate = cmd_ln_int32_r(bd->config, "-frate");
    int16 buf[512];

    if (ef != -1) {
        ef = (int32) (((ef - sf) * samprate / frate) + (samprate
                * cmd_ln_float32_r(bd->config, "-wlen")));
    }
    sf = (int32) (sf * (samprate / frate));
    fseek(infh, cmd_ln_int32_r(bd->config, "-adchdr") + sf * sizeof(int16),
            SEEK_SET);

    if (al) {
        alignment_iter_t *itor;
        double starttime = 0.0;
        for (itor = alignment_words(al); itor; itor = alignment_iter_next(itor)) {
            alignment_entry_t *ent = alignment_iter_get(itor);
            double nsec = (double) ent->duration / frate;
            double endtime = starttime + nsec;
            size_t nsamp = (size_t) (nsec * samprate);
            E_INFO("Processing %d samples for %s (%f seconds ending %f)\n",
            nsamp, dict_wordstr(search_factory_d2p(bd->sf)->dict, ent->id.wid),
            nsec, endtime);
            E_INFO("Woke up at delta %f\n", get_time_delta(bd));
            while (nsamp > 0) {
                size_t nread = 512;
                if (nread > nsamp) nread = nsamp;
                nread = fread(buf, sizeof(int16), nread, infh);
                if (nread == 0)
                break;
                featbuf_producer_process_raw(fb, buf, nread, FALSE);
                nsamp -= nread;
                starttime += (nread / samprate);
                double delta = get_time_delta(bd);
                if (starttime > delta) {
                    E_INFO("Sleeping until next start time (%f seconds)\n",
                            starttime - delta);
                    usleep((int)((starttime - delta) * 1000000));
                }
            }
            double delta = get_time_delta(bd);
            if (endtime > delta) {
                E_INFO("Sleeping until end time (%f seconds)\n",
                        endtime - delta);
                usleep((int)((endtime - delta) * 1000000));
            }
        }
    }
    else {
        while (ef == -1 || sf < ef) {
            size_t nread = 512;
            if (ef != -1 && nread > ef - sf)
            nread = ef - sf;
            nread = fread(buf, sizeof(int16), nread, infh);
            if (nread == 0)
            break;
            featbuf_producer_process_raw(fb, buf, nread, FALSE);
            //usleep((int)((double)nread / 16000 * 1000000));
            sf += nread;
        }
    }
    return 0;
}

static int batch_decoder_decode_mfc(batch_decoder_t *bd, FILE *infh, int sf,
        int ef, alignment_t *al)
{
    featbuf_t *fb = search_factory_featbuf(bd->sf);
    mfcc_t **mfcs;
    int nfr, rv;

    if (NULL == (mfcs = read_mfc_file(infh, sf, ef, &nfr,
            cmd_ln_int32_r(bd->config, "-ceplen"))))
        return -1;

    rv = featbuf_producer_process_cep(fb, mfcs, nfr, TRUE);
    ckd_free_2d(mfcs);
    return rv;
}

int batch_decoder_decode(batch_decoder_t *bd, char *file, char *uttid,
        int32 sf, int32 ef, alignment_t *al)
{
    featbuf_t *fb;
    FILE *infh;
    char const *cepdir, *cepext;
    char *infile;
    int rv;

    if (ef != -1 && ef < sf) {
        E_ERROR("End frame %d is < start frame %d\n", ef, sf);
        return -1;
    }

    cepdir = cmd_ln_str_r(bd->config, "-cepdir");
    cepext = cmd_ln_str_r(bd->config, "-cepext");

    /* Build input filename. */
    infile = string_join(cepdir ? cepdir : "",
            "/", file,
            cepext ? cepext : "", NULL);
    if (uttid == NULL) uttid = file;

    if ((infh = fopen(infile, "rb")) == NULL)
    {
        E_ERROR_SYSTEM("Failed to open %s", infile);
        return -1;
    }

    fb = search_factory_featbuf(bd->sf);
    gettimeofday(&bd->utt_start, NULL);
    featbuf_producer_start_utt(fb, uttid);

    if (cmd_ln_boolean_r(bd->config, "-adcin"))
    rv = batch_decoder_decode_adc(bd, infh, sf, ef, al);
    else
    rv = batch_decoder_decode_mfc(bd, infh, sf, ef, al);

    featbuf_producer_end_utt(fb);
    if (bd->hypfh) {
        char const *hyp;
        int32 score;
        hyp = search_hyp(bd->fwdflat, &score);
        fprintf(bd->hypfh, "%s (%s %d)\n",
                hyp, uttid, score);
    }

    fclose(infh);
    ckd_free(infile);

    return rv;
}

alignment_t *
parse_alignment(char *line, dict2pid_t *d2p)
{
    alignment_t *al;
    char **wptr;
    int nf, i;
    double spos;
    int32 frate = 100; /* FIXME */

    nf = str2words(line, NULL, 0);
    if (nf < 0)
        return NULL;
    wptr = ckd_calloc(nf, sizeof(*wptr));
    nf = str2words(line, wptr, nf);
    if (nf < 0) {
        ckd_free(wptr);
        return NULL;
    }
    al = alignment_init(d2p);
    spos = 0.0;
    for (i = 0; i < nf; ++i) {
        char *c = strchr(wptr[i], ':');
        double epos;
        int duration;
        if (c == NULL) /* word ID */
            break;
        *c++ = '\0';
        epos = atof(c);
        duration = (int) ((epos - spos) * frate);
        alignment_add_word(al, dict_wordid(d2p->dict, wptr[i]), duration);
        spos = epos;
    }
    return al;
}

int batch_decoder_run(batch_decoder_t *bd)
{
    int32 ctloffset, ctlcount, ctlincr;
    lineiter_t *li, *ali = NULL;

    search_run(bd->fwdtree);
    search_run(bd->fwdflat);

    ctloffset = cmd_ln_int32_r(bd->config, "-ctloffset");
    ctlcount = cmd_ln_int32_r(bd->config, "-ctlcount");
    ctlincr = cmd_ln_int32_r(bd->config, "-ctlincr");

    if (bd->alignfh)
        ali = lineiter_start(bd->alignfh);
    for (li = lineiter_start(bd->ctlfh); li; li = lineiter_next(li)) {
        alignment_t *al = NULL;
        char *wptr[4];
        int32 nf, sf, ef;

        if (li->lineno < ctloffset) {
            if (ali)
                ali = lineiter_next(ali);
            continue;
        }
        if ((li->lineno - ctloffset) % ctlincr != 0) {
            if (ali)
                ali = lineiter_next(ali);
            continue;
        }
        if (ctlcount != -1 && li->lineno >= ctloffset + ctlcount)
            break;
        if (ali)
            al = parse_alignment(ali->buf, search_factory_d2p(bd->sf));
        sf = 0;
        ef = -1;
        nf = str2words(li->buf, wptr, 4);
        if (nf == 0) {
            /* Do nothing. */
        }
        else if (nf < 0) {
            E_ERROR("Unexpected extra data in control file at line %d\n", li->lineno);
        }
        else
        {
            char *file, *uttid;
            file = wptr[0];
            uttid = NULL;
            if (nf > 1)
            sf = atoi(wptr[1]);
            if (nf > 2)
            ef = atoi(wptr[2]);
            if (nf > 3)
            uttid = wptr[3];
            /* Do actual decoding. */
            batch_decoder_decode(bd, file, uttid, sf, ef, al);
        }
        alignment_free(al);
        if (ali) ali = lineiter_next(ali);
    }
    featbuf_producer_shutdown(search_factory_featbuf(bd->sf));
    return 0;
}

static int search_cb(search_t *search, search_event_t *evt, void *udata)
{
    batch_decoder_t *bd = (batch_decoder_t *) udata;
    dict_t *d = search_factory_d2p(bd->sf)->dict;
    double delta = get_time_delta(bd);
    double frate = cmd_ln_int32_r(search_config(search), "-frate");
    FILE *hypfh = NULL;
    void *val;

    if (hash_table_lookup(bd->hypfiles, search_name(search), &val) == 0)
        hypfh = val;
    else
        hypfh = stdout;
    fprintf(hypfh, "time delta %f ", delta);
    switch (evt->event) {
        case SEARCH_PARTIAL_RESULT: {
            int32 score;
            seg_iter_t *seg = search_seg_iter(search, &score);
            fprintf(hypfh, "partial: ");
            for (; seg; seg = seg_iter_next(seg)) {
                int sf, ef;
                seg_iter_times(seg, &sf, &ef);
                fprintf(hypfh, "%s:%.3f ", dict_basestr(d, seg_iter_wid(seg)),
                        (double) ef / frate);
            }
            fprintf(hypfh, "(%s)\n", search_uttid(search));
            break;
        }
        case SEARCH_START_UTT:
            fprintf(hypfh, "start %s\n", search_uttid(search));
            break;
        case SEARCH_END_UTT:
            fprintf(hypfh, "end %s\n", search_uttid(search));
            break;
        case SEARCH_FINAL_RESULT: {
            int32 score;
            seg_iter_t *seg = search_seg_iter(search, &score);
            fprintf(hypfh, "full: ");
            for (; seg; seg = seg_iter_next(seg)) {
                int sf, ef;
                seg_iter_times(seg, &sf, &ef);
                fprintf(hypfh, "%s:%.3f ", dict_basestr(d, seg_iter_wid(seg)),
                        (double) ef / frate);
            }
            fprintf(hypfh, "(%s)\n", search_uttid(search));
            break;
        }
    }
    return 0;
}

batch_decoder_t *
batch_decoder_init_argv(int argc, char *argv[])
{
    batch_decoder_t *bd;
    char const *str;

    bd = ckd_calloc(1, sizeof(*bd));
    bd->config = cmd_ln_parse_r(NULL, ms_args_def, argc, argv, FALSE);
    if ((str = cmd_ln_str_r(bd->config, "-ctl")) == NULL) {
        E_ERROR("-ctl argument not present, nothing to do in batch mode!\n");
        goto error_out;
    }
    if ((bd->ctlfh = fopen(str, "r")) == NULL) {
        E_ERROR_SYSTEM("Failed to open control file '%s'", str);
        goto error_out;
    }
    if ((str = cmd_ln_str_r(bd->config, "-align")) != NULL) {
        if ((bd->alignfh = fopen(str, "r")) == NULL) {
            E_ERROR_SYSTEM("Failed to open align file '%s'", str);
        }
    }
    if ((str = cmd_ln_str_r(bd->config, "-hyp")) != NULL) {
        if ((bd->hypfh = fopen(str, "w")) == NULL) {
            E_ERROR_SYSTEM("Failed to open hypothesis file '%s'", str);
        }
    }

    if ((bd->sf = search_factory_init_cmdln(bd->config)) == NULL)
        goto error_out;
    if ((str = cmd_ln_str_r(bd->config, "-fwdtreelm")) != NULL) {
        if ((bd->fwdtree = search_factory_create(bd->sf, NULL, "fwdtree",
                "-fwdtreelm", str, NULL)) == NULL)
            goto error_out;
        if ((bd->fwdflat = search_factory_create(bd->sf, NULL, "fwdflat", NULL)) == NULL)
            goto error_out;
    }
    else {
        if ((bd->fwdtree = search_factory_create(bd->sf, NULL, "fwdtree", NULL)) == NULL)
            goto error_out;
        if ((bd->fwdflat = search_factory_create(bd->sf, bd->fwdtree, "fwdflat", NULL)) == NULL)
            goto error_out;
    }
    if ((str = cmd_ln_str_r(bd->config, "-vm")) != NULL) {
        vocab_map_t *vm = vocab_map_init(search_factory_d2p(bd->sf)->dict);
        FILE *vmfh;
        if (vm == NULL)
            goto error_out;
        if ((vmfh = fopen(str, "r")) == NULL) {
            vocab_map_free(vm);
            goto error_out;
        }
        if (vocab_map_read(vm, vmfh) < 0) {
            vocab_map_free(vm);
            goto error_out;
        }
        fclose(vmfh);
        fwdflat_search_set_vocab_map(bd->fwdflat, vm);
    }
    //if ((bd->latgen = search_factory_create(bd->sf, "latgen", NULL)) == NULL)
    //goto error_out;

    search_link(bd->fwdtree, bd->fwdflat, "fwdtree", FALSE);
    // search_link(bd->fwdflat, bd->latgen, "fwdflat", TRUE);
    search_set_cb(bd->fwdtree, search_cb, bd);
    search_set_cb(bd->fwdflat, search_cb, bd);

    bd->hypfiles = hash_table_new(0, FALSE);
    if ((str = cmd_ln_str_r(bd->config, "-hypprefix"))) {
        char *hypfile;
        FILE *hypfh;
        hypfile = string_join(str, ".fwdtree.hyp", NULL);
        hypfh = fopen(hypfile, "w");
        if (hypfh == NULL) {
            E_ERROR_SYSTEM("Could not open %s", hypfile);
        }
        else {
            hash_table_enter(bd->hypfiles, "fwdtree", hypfh);
        }
        ckd_free(hypfile);
        hypfile = string_join(str, ".fwdflat.hyp", NULL);
        hypfh = fopen(hypfile, "w");
        if (hypfh == NULL) {
            E_ERROR_SYSTEM("Could not open %s", hypfile);
        }
        else {
            hash_table_enter(bd->hypfiles, "fwdflat", hypfh);
        }
        ckd_free(hypfile);
    }
    return bd;

    error_out:
    return NULL;
}

int batch_decoder_free(batch_decoder_t *bd)
{
    hash_iter_t *itor;
    if (bd == NULL)
        return 0;
    if (bd->ctlfh != NULL)
        fclose(bd->ctlfh);
    if (bd->alignfh != NULL)
        fclose(bd->alignfh);
    if (bd->hypfh != NULL)
        fclose(bd->hypfh);
    cmd_ln_free_r(bd->config);
    search_free(bd->fwdtree);
    search_free(bd->fwdflat);
    //search_free(bd->latgen);
    search_factory_free(bd->sf);
    for (itor = hash_table_iter(bd->hypfiles); itor; itor
            = hash_table_iter_next(itor)) {
        fclose(hash_entry_val(itor->ent));
    }
    hash_table_free(bd->hypfiles);
    ckd_free(bd);
    return 0;
}

int main(int argc, char *argv[])
{
    batch_decoder_t *bd;

    if ((bd = batch_decoder_init_argv(argc, argv)) == NULL) {
        E_ERROR("Failed to initialize decoder\n");
        return 1;
    }

    batch_decoder_run(bd);
    batch_decoder_free(bd);

    return 0;
}
