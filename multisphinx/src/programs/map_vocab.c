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

/* PocketSphinx headers. */
#include <pocketsphinx.h>

static const arg_t ps_args_def[] = {
    POCKETSPHINX_OPTIONS,
    { "-lsn",
      ARG_STRING,
      NULL,
      "Input transcription file." },
    { "-bgdict",
      ARG_STRING,
      NULL,
      "Background dictionary file." },
    { "-prune_topn",
      ARG_INT,
      NULL,
      "In absence of -bgdict, prune vocabulary by taking top-N unigrams." },
    /* Other methods will certainly folllow */
    CMDLN_EMPTY_OPTION
};

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
						   -s->mixw[j][k][senid]);
			for (i = 0; i < s->n_sen; ++i) {
				acmod->senone_scores[i]
					+= (int)(mixw * (-s->mixw[j][k][senid]
							 - -s->mixw[j][k][i]))
					/* fake_kl_scaled */
					/ s->n_feat;
			}
		}
	}
}

int
main(int argc, char *argv[])
{
	return 0;
}
