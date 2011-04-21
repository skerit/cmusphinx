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
 * @file ps_alignment.h Multi-level alignment structure
 */

#ifndef __PS_ALIGNMENT_H__
#define __PS_ALIGNMENT_H__

/* System headers. */
#include <stdarg.h>

/* SphinxBase headers. */
#include <sphinxbase/prim_type.h>

/* MultiSphinx headers. */
#include <multisphinx/dict2pid.h>
#include <multisphinx/hmm.h>

#define PS_ALIGNMENT_NONE ((uint16)0xffff)

struct alignment_entry_s {
    union {
        int32 wid;
        struct {
            uint16 ssid;
            uint16 cipid;
            uint16 tmatid;
        } pid;
        uint16 senid;
    } id;
    int16 start;
    int16 duration;
    uint16 parent;
    uint16 child;
};
typedef struct alignment_entry_s alignment_entry_t;

struct alignment_vector_s {
    alignment_entry_t *seq;
    uint16 n_ent, n_alloc;
};
typedef struct alignment_vector_s alignment_vector_t;

struct alignment_s {
    int refcnt;
    dict2pid_t *d2p;
    alignment_vector_t word;
    alignment_vector_t sseq;
    alignment_vector_t state;
};
typedef struct alignment_s alignment_t;

struct alignment_iter_s {
    alignment_t *al;
    alignment_vector_t *vec;
    int pos;
};
typedef struct alignment_iter_s alignment_iter_t;

/**
 * Create a new, empty alignment.
 */
alignment_t *alignment_init(dict2pid_t *d2p);

/**
 * Retain an alignment
 */
alignment_t *alignment_retain(alignment_t *al);

/**
 * Release an alignment
 */
int alignment_free(alignment_t *al);

/**
 * Append a word to the alignment.
 *
 * @param al Alignment to add the word to.
 * @param wid Word ID for the word.
 * @param duration Duration of the word (or zero if no explicit duration given).
 */
int alignment_add_word(alignment_t *al,
                       int32 wid, int duration);

/**
 * Append a sequence of words to the alignment using its internal dictionary.
 */
int alignment_add_words(alignment_t *al, char const *w1, ...);

/**
 * Populate lower layers using available word information.
 */
int alignment_populate(alignment_t *al);

/**
 * Populate lower layers using context-independent phones.
 */
int alignment_populate_ci(alignment_t *al);

/**
 * Propagate timing information up from state sequence.
 */
int alignment_propagate(alignment_t *al);

/**
 * Number of words.
 */
int alignment_n_words(alignment_t *al);

/**
 * Number of phones.
 */
int alignment_n_phones(alignment_t *al);

/**
 * Number of states.
 */
int alignment_n_states(alignment_t *al);

/**
 * Iterate over the alignment starting at the first word.
 */
alignment_iter_t *alignment_words(alignment_t *al);

/**
 * Iterate over the alignment starting at the first phone.
 */
alignment_iter_t *alignment_phones(alignment_t *al);

/**
 * Iterate over the alignment starting at the first state.
 */
alignment_iter_t *alignment_states(alignment_t *al);

/**
 * Get the alignment entry pointed to by an iterator.
 */
alignment_entry_t *alignment_iter_get(alignment_iter_t *itor);

/**
 * Move alignment iterator to given index.
 */
alignment_iter_t *alignment_iter_goto(alignment_iter_t *itor, int pos);

/**
 * Move an alignment iterator forward.
 */
alignment_iter_t *alignment_iter_next(alignment_iter_t *itor);

/**
 * Move an alignment iterator back.
 */
alignment_iter_t *alignment_iter_prev(alignment_iter_t *itor);

/**
 * Get a new iterator starting at the parent of the current node.
 */
alignment_iter_t *alignment_iter_up(alignment_iter_t *itor);
/**
 * Get a new iterator starting at the first child of the current node.
 */
alignment_iter_t *alignment_iter_down(alignment_iter_t *itor);

/**
 * Release an iterator before completing all iterations.
 */
int alignment_iter_free(alignment_iter_t *itor);

#endif /* __PS_ALIGNMENT_H__ */
