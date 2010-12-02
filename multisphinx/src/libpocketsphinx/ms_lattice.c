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
 * @file ms_lattice.c Word lattices for MultiSphinx.
 */

#include <sphinxbase/pio.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/ckd_alloc.h>

#include "ms_lattice.h"

ms_lattice_t *
ms_lattice_init(logmath_t *lmath)
{
    ms_lattice_t *l = ckd_calloc(1, sizeof(*l));
    l->refcount = 1;
    l->lmath = logmath_retain(lmath);
    l->node_list = garray_init(0, sizeof(ms_latnode_t));
    l->node_map = nodeid_map_init();
    return l;
}

ms_lattice_t *
ms_lattice_retain(ms_lattice_t *l)
{
    ++l->refcount;
    return l;
}

int
ms_lattice_free(ms_lattice_t *l)
{
    if (l == NULL)
        return 0;
    if (--l->refcount > 0)
        return l->refcount;

    garray_free(l->node_list);
    nodeid_map_free(l->node_map);
    logmath_free(l->lmath);
    ckd_free(l);
    return 0;
}

int32
ms_lattice_node(ms_lattice_t *l, int32 lmstate, int sf)
{
}

ms_latnode_t *
ms_lattice_get_node_idx(ms_lattice_t *l, int32 idx)
{
}

ms_latnode_t *
ms_lattice_get_node_id(ms_lattice_t *l, nodeid_t id)
{
}

ms_latnode_t *
ms_lattice_get_start(ms_lattice_t *l)
{
}

ms_latnode_t *
ms_lattice_get_end(ms_lattice_t *l)
{
}

ms_latlink_t *
ms_lattice_link(ms_lattice_t *l,
                ms_latnode_t *src, ms_latnode_t *dest,
                int32 wid, int32 ascr)
{
}

int
ms_lattice_read_htk(ms_lattice_t *l, FILE *fh)
{
}

int
ms_lattice_write_htk(ms_lattice_t *l, FILE *fh)
{
}

int
ms_lattice_write_dot(ms_lattice_t *l, FILE *fh)
{
}
