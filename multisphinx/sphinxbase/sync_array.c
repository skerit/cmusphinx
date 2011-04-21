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
 * @file sync_array.c Expandable arrays with synchronization.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include <time.h>
#include <string.h>

#include <sphinxbase/sbthread.h>
#include <sphinxbase/garray.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/err.h>

#include "sync_array.h"

struct sync_array_s {
    int refcount;
    garray_t *data;
    garray_t *count;
    sbmtx_t *mtx;
    sbevent_t *evt;
    size_t final_next_idx;
};

sync_array_t *
sync_array_init(size_t n_ent, size_t ent_size)
{
    sync_array_t *sa;

    sa = ckd_calloc(1, sizeof(*sa));
    sa->refcount = 1;
    sa->data = garray_init(n_ent, ent_size);
    sa->count = garray_init(n_ent, 1);
    sa->mtx = sbmtx_init();
    sa->evt = sbevent_init(FALSE);
    sa->final_next_idx = (size_t)-1;

    return sa;
}

sync_array_t *
sync_array_retain(sync_array_t *sa)
{
    sbmtx_lock(sa->mtx);
    if (sa->refcount == 255) {
        E_INFO("Failed to retain sync_array %p, refcount has reached 255\n", sa);
        sbmtx_unlock(sa->mtx);
        return NULL;
    }
    ++sa->refcount;
    sbmtx_unlock(sa->mtx);
    return sa;
}

int
sync_array_free(sync_array_t *sa)
{
    if (sa == NULL)
        return 0;
    sbmtx_lock(sa->mtx);
    if (--sa->refcount > 0) {
        int refcount = sa->refcount;
        size_t i;
        /* FIXME: This may lead to memory leaks.  We don't know
         * exactly which elements this thread had laid claim to, so we
         * have to decrement the count on all of them.  Best practice,
         * as described in the header, is to have a consumer release
         * all remaining elements before freeing the array. */
        for (i = garray_base(sa->count);
             i < garray_next_idx(sa->count); ++i)
            if (garray_ent(sa->count, uint8, i) > 0)
                --garray_ent(sa->count, uint8, i);
        sbmtx_unlock(sa->mtx);
        return refcount;
    }		
    sbmtx_unlock(sa->mtx);
    E_INFO("Maximum allocation %d items (%d KiB)\n",
           garray_alloc_size(sa->data),
           garray_alloc_size(sa->data) * (garray_ent_size(sa->data) + 1) / 1024);
    garray_free(sa->data);
    garray_free(sa->count);
    sbevent_free(sa->evt);
    sbmtx_free(sa->mtx);
    ckd_free(sa);
    return 0;
}

size_t
sync_array_next_idx(sync_array_t *sa)
{
    return garray_next_idx(sa->data);
}

int
sync_array_wait(sync_array_t *sa, size_t idx, int sec, int nsec)
{
    int tsec, tnsec, nwait = 0;

    if (sec == -1) {
        tsec = 0;
        tnsec = 50000; /* Arbitrary polling interval. */
    }
    else {
        tsec = sec;
        tnsec = nsec;
    }

    /* Wait until length of sa->data > idx or end of utt. */
    while (1) {
        sbmtx_lock(sa->mtx);
        if (garray_next_idx(sa->data) > idx) {
            sbmtx_unlock(sa->mtx);
            return 0;
        }
        if (idx >= sa->final_next_idx) {
            E_INFO("idx %d is final (%d)\n", idx, sa->final_next_idx);
            sbmtx_unlock(sa->mtx);
            return -1;
        }
        sbmtx_unlock(sa->mtx);
        if (nwait > 0)
            return -1;
        /* If we wait forever here, there's a race condition
         * between the tests above and the wait here.  So we
         * poll the array no matter what, and if we had a
         * timeout, we make sure to only do it once. */
        if (sbevent_wait(sa->evt, tsec, tnsec) < 0)
            return -1;
        if (sec != -1)
            ++nwait;
    }
    /* Never reached. */
    return -1;
}

int
sync_array_get(sync_array_t *sa, size_t idx, void *out_ent)
{
    sbmtx_lock(sa->mtx);
    if (idx < garray_base(sa->data)
        || idx >= garray_next_idx(sa->data)) {
        sbmtx_unlock(sa->mtx);
        return -1;
    }
    memcpy(out_ent, garray_void(sa->data, idx),
           garray_ent_size(sa->data));
    sbmtx_unlock(sa->mtx);
    return 0;
}


int
sync_array_append(sync_array_t *sa, void *ent)
{
    int zero = 0;
    void *new_ent;

    sbmtx_lock(sa->mtx);
    /* Not allowed to append to a finalized array. */
    if (garray_next_idx(sa->data) >= sa->final_next_idx) {
        sbmtx_unlock(sa->mtx);
        return -1;
    }
    new_ent = garray_append(sa->data, ent);
    garray_append(sa->count, &zero);
    sbevent_signal(sa->evt);
    sbmtx_unlock(sa->mtx);

    return 0;
}

size_t
sync_array_finalize(sync_array_t *sa)
{
    sbmtx_lock(sa->mtx);
    /* Not allowed to do this more than once! (or from multiple
     * threads at the same time) */
    if (sa->final_next_idx != (size_t) -1) {
        sbmtx_unlock(sa->mtx);
        return -1;
    }
    sa->final_next_idx = garray_next_idx(sa->data);
    sbmtx_unlock(sa->mtx);

    return sa->final_next_idx;
}

int
sync_array_force_quit(sync_array_t *sa)
{
    sbmtx_lock(sa->mtx);
    /* Guaranteed to make everything fail. */
    sa->final_next_idx = 0;
    sbmtx_unlock(sa->mtx);
    return 0;
}

int
sync_array_reset(sync_array_t *sa)
{
    sbmtx_lock(sa->mtx);
    garray_reset(sa->data);
    garray_reset(sa->count);
    sa->final_next_idx = (size_t)-1;
    sbmtx_unlock(sa->mtx);
    return 0;
}

size_t
sync_array_release(sync_array_t *sa, size_t start_idx, size_t end_idx)
{
    size_t i;

    sbmtx_lock(sa->mtx);
    if (start_idx < garray_base(sa->count))
        start_idx = garray_base(sa->count);
    if (start_idx > garray_next_idx(sa->count))
        start_idx = garray_next_idx(sa->count);
    if (end_idx > garray_next_idx(sa->count))
        end_idx = garray_next_idx(sa->count);
    if (end_idx <= start_idx) {
        sbmtx_unlock(sa->mtx);
        return start_idx;
    }
    /* Increment count for all indices. */
    for (i = start_idx; i < end_idx; ++i)
        ++garray_ent(sa->count, uint8, i);

    /* Print stuff for debugging. */
#if 0
    printf("rc %d counts[%d:%d]:", (int)sa->refcount,
           (int)garray_base(sa->count),
           (int)garray_next_idx(sa->count));
    for (i = garray_base(sa->count);
         i < garray_next_idx(sa->count); ++i)
        printf(" %d", (int)garray_ent(sa->count, uint8, i));
    printf("\n");
#endif

    /* Find first reachable element. */
    for (i = garray_base(sa->count);
         i < garray_next_idx(sa->count); ++i)
        /* Note that we assume the producer retains one
         * reference to the array. */
        if (garray_ent(sa->count, uint8, i) < sa->refcount - 1)
            break;

    /* Release unreachable elements. */
    if (i > garray_base(sa->count)) {
        /* printf("Releasing up to %d\n", (int)i); */
        garray_shift_from(sa->count, i);
        garray_shift_from(sa->data, i);
        garray_set_base(sa->count, i);
        garray_set_base(sa->data, i);
    }
    sbmtx_unlock(sa->mtx);

    return i;
}

size_t
sync_array_release_all(sync_array_t *sa)
{
    return sync_array_release(sa, 0, (size_t)-1);
}

size_t
sync_array_available(sync_array_t *sa)
{
    return garray_base(sa->count);
}
