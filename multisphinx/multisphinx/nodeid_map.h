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

#ifndef __NODEID_MAP_H__
#define __NODEID_MAP_H__

/**
 * @file nodeid_map.h Node ID maps.
 */

#include <sphinxbase/prim_type.h>

/**
 * Lattice node ID.
 */
typedef struct nodeid_s {
    int32 lmstate; /**< Language model state ID. */
    int16 sf;   /**< Start frame. */
} nodeid_t;

/**
 * Map of node IDs (LM state/frame pairs) to lattice indices.
 */
typedef struct nodeid_map_s nodeid_map_t;

/**
 * Iterator over node IDs
 */
typedef struct nodeid_iter_s nodeid_iter_t;

nodeid_map_t *nodeid_map_init(void);
int nodeid_map_free(nodeid_map_t *nmap);
int32 nodeid_map_add(nodeid_map_t *nmap, int sf, int32 lmstate, int32 idx);
int32 nodeid_map_remap(nodeid_map_t *nmap, int sf, int32 lmstate, int32 idx);
int32 nodeid_map_delete(nodeid_map_t *nmap, int sf, int32 lmstate);
int32 nodeid_map_map(nodeid_map_t *nmap, int sf, int32 lmstate);

/**
 * Iterate over node IDs in frame order.
 *
 * If no start frame is given then all nodes will be iterated over in
 * frame order.
 */
nodeid_iter_t *nodeid_map_iter(nodeid_map_t *nmap, int sf);
nodeid_iter_t *nodeid_iter_next(nodeid_iter_t *itor);
int32 nodeid_iter_get(nodeid_iter_t *itor, int32 *out_lmstate);
void nodeid_iter_free(nodeid_iter_t *itor);

#endif /* __NODEID_MAP_H__ */
