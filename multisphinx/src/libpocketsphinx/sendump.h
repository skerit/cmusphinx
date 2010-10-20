/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 1999-2010 Carnegie Mellon University.  All rights
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
 * @file sendump.h Senone dump file
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __SENDUMP_H__
#define __SENDUMP_H__

typedef struct sendump_s {
    int refcount;
    uint8 *sen2cb;     /**< Senone to codebook mapping. */
    uint8 ***mixw;     /**< Mixture weight distributions by feature, codeword, senone */
    mmio_file_t *sendump_mmap;/* Memory map for mixw (or NULL if not mmap) */
    uint8 *mixw_cb;    /* Mixture weight codebook, if any (assume it contains 16 values) */
} sendump_t;

sendump_t *sendump_read_sendump(cmd_ln_t *config, logmath_t *lmath_8b,
                                gauden_t *g, bin_mdef_t *mdef,
                                char const *file);


sendump_t *sendump_read_mixw(cmd_ln_t *config, logmath_t *lmath_8b,
                             gauden_t *g, bin_mdef_t *mdef,
                             char const *file);

sendump_t *sendump_retain(sendump_t *s);

int sendump_free(sendump_t *s);

#endif /* __SENDUMP_H__ */
