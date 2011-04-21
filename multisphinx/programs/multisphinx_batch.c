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

typedef struct batch_decoder_s
{
    search_factory_t *sf;
    cmd_ln_t *config;
    search_t *fwdtree;
    search_t *fwdflat;
    search_t *latgen;

    FILE *ctlfh;
} batch_decoder_t;

/**
 * Command-line argument definitions for batch processing.
 */
static const arg_t ms_args_def[] =
{ POCKETSPHINX_OPTIONS,
/* Various options specific to batch-mode processing. */
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

static int batch_decoder_decode_adc(batch_decoder_t *bd, FILE *infh, int sf, int ef)
{
    featbuf_t *fb = search_factory_featbuf(bd->sf);
    int16 *buf;

    if (ef != -1)
    {
        ef = (int32) ((ef - sf) * (cmd_ln_float32_r(bd->config, "-samprate")
                / cmd_ln_int32_r(bd->config, "-frate"))
                + (cmd_ln_float32_r(bd->config, "-samprate")
                        * cmd_ln_float32_r(bd->config, "-wlen")));
    }
    sf = (int32) (sf * (cmd_ln_float32_r(bd->config, "-samprate")
            / cmd_ln_int32_r(bd->config, "-frate")));
    fseek(infh, cmd_ln_int32_r(bd->config, "-adchdr") + sf * sizeof(int16),
            SEEK_SET);

    // FIXME: should try to do whole utterance if possible...
    while (sf < ef) {
        size_t nread = fread(buf, sizeof(int16), 2048, infh);
        if (nread == 0)
            break;
        if (nread > ef - sf)
            nread = ef - sf;
        featbuf_producer_process_raw(fb, buf, nread, FALSE);
        sf += nread;
    }
    return 0;
}

static int batch_decoder_decode_mfc(batch_decoder_t *bd, FILE *infh, int sf, int ef)
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

int batch_decoder_decode(batch_decoder_t *bd, char *file,
        char *uttid, int32 sf, int32 ef)
{
    featbuf_t *fb;
    FILE *infh;
    char const *cepdir, *cepext;
    char *infile;
    int rv;

    if (ef != -1 && ef < sf)
    {
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
        ckd_free(infile);
        return -1;
    }

    fb = search_factory_featbuf(bd->sf);
    featbuf_producer_start_utt(fb, uttid);

    if (cmd_ln_boolean_r(bd->config, "-adcin"))
        rv = batch_decoder_decode_adc(bd, infh, sf, ef);
    else
        rv = batch_decoder_decode_mfc(bd, infh, sf, ef);

    /* FIXME: This will wait for the utterance to be over, which isn't actually what we want. */
    featbuf_producer_end_utt(fb);

    fclose(infh);
    ckd_free(infile);

    return rv;
}

int batch_decoder_run(batch_decoder_t *bd)
{
    int32 ctloffset, ctlcount, ctlincr;
    lineiter_t *li;

    search_run(bd->fwdtree);
    search_run(bd->fwdflat);

    ctloffset = cmd_ln_int32_r(bd->config, "-ctloffset");
    ctlcount = cmd_ln_int32_r(bd->config, "-ctlcount");
    ctlincr = cmd_ln_int32_r(bd->config, "-ctlincr");

    for (li = lineiter_start(bd->ctlfh); li; li = lineiter_next(li))
    {
        char *wptr[4];
        int32 nf, sf, ef;

        if (li->lineno < ctloffset)
            continue;
        if ((li->lineno - ctloffset) % ctlincr != 0)
            continue;
        if (ctlcount != -1 && li->lineno >= ctloffset + ctlcount)
            break;
        sf = 0;
        ef = -1;
        nf = str2words(li->buf, wptr, 4);
        if (nf == 0)
        {
            /* Do nothing. */
        }
        else if (nf < 0)
        {
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
            batch_decoder_decode(bd, file, uttid, sf, ef);
        }
    }
    return 0;
}

batch_decoder_t *
batch_decoder_init_argv(int argc, char *argv[])
{
    batch_decoder_t *bd;
    char const *ctl;

    bd = ckd_calloc(1, sizeof(*bd));
    bd->config = cmd_ln_parse_r(NULL, ms_args_def, argc, argv, FALSE);
    if ((ctl = cmd_ln_str_r(bd->config, "-ctl")) == NULL)
    {
        E_ERROR("-ctl argument not present, nothing to do in batch mode!\n");
        goto error_out;
    }
    if ((bd->ctlfh = fopen(ctl, "r")) == NULL)
    {
        E_ERROR_SYSTEM("Failed to open control file '%s'", ctl);
        goto error_out;
    }

    if ((bd->sf = search_factory_init_cmdln(bd->config)) == NULL)
    goto error_out;
    if ((bd->fwdtree = search_factory_create(bd->sf, "fwdtree", NULL)) == NULL)
    goto error_out;
    if ((bd->fwdflat = search_factory_create(bd->sf, "fwdflat", NULL)) == NULL)
    goto error_out;
    if ((bd->latgen = search_factory_create(bd->sf, "latgen", NULL)) == NULL)
    goto error_out;

    search_link(bd->fwdtree, bd->fwdflat, "fwdtree", FALSE);
    // search_link(bd->fwdflat, bd->latgen, "fwdflat", TRUE);

    return bd;

    error_out:
    return NULL;
}

int batch_decoder_free(batch_decoder_t *bd)
{
    if (bd == NULL)
        return 0;
    if (bd->ctlfh != NULL)
        fclose(bd->ctlfh);
    cmd_ln_free_r(bd->config);
    search_free(bd->fwdtree);
    search_free(bd->fwdflat);
    search_free(bd->latgen);
    search_factory_free(bd->sf);
    ckd_free(bd);
    return 0;
}

int main(int argc, char *argv[])
{
    batch_decoder_t *bd;

    if ((bd = batch_decoder_init_argv(argc, argv)) == NULL)
    {
        E_ERROR("Failed to initialize decoder");
        return 1;
    }

    batch_decoder_run(bd);

    batch_decoder_free(bd);

    return 0;
}
