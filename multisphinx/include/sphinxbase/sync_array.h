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
 * @file sync_array.h Expandable arrays with synchronization.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __SYNC_ARRAY_H__
#define __SYNC_ARRAY_H__

#include <stdlib.h>

/**
 * Expandable, bounds-checked array with synchronization.
 *
 * This is an implementation of a growable array with a movable base
 * pointer, similar to garray_t, with the added functionality of
 * synchronization between a single producer and multiple consumers.
 * The producer can append to the end of the array, while the
 * consumers have random access to its contents.
 *
 * In addition, the elements of the array are reference counted and
 * can be released by consumers.  When all consumers have released
 * claims on an initial sequence of the array, the memory associated
 * with it will be released.  Since this implies that the elements may
 * be moved in memory, the array cannot be accessed through pointers.
 */
typedef struct sync_array_s sync_array_t;

/**
 * Create and initialize a new sync_array_t.
 *
 * @param n_ent Initial number of elements.
 * @param ent_size Size of each element.
 * @return Newly initialized array, or NULL on failure.
 */
sync_array_t *sync_array_init(size_t n_ent, size_t ent_size);

/**
 * Retain a pointer to a sync_array_t.
 *
 * In order for element reference counting and release to work
 * properly this function must be called once for each consumer
 * thread.
 *
 * It is possible for this to fail. In particular, only 254 consumer
 * threads are allowed.  In this case it returns NULL.
 *
 * @param sa Array.
 * @return Pointer to array, with reference count increased, or NULL
 *         on failure.
 */
sync_array_t *sync_array_retain(sync_array_t *sa);

/**
 * Release a pointer to a sync_array_t.
 *
 * Note that this DOES NOT release any elements of the array that a
 * consumer might have claim to.  Best practice is for consumer
 * threads to release all remaining elements of the array before
 * calling this, in order to preserve the memory efficiency of the
 * array.
 *
 * @param sa Array.
 * @return New reference count.
 */
int sync_array_free(sync_array_t *sa);

/**
 * Get the index of the first available element in an array.
 *
 * This function returns the index of the first element which has not
 * yet been released by all consumers.  It can also be used to detect
 * when all consumers have finished, as the index returned will then
 * be equal to the one returned by sync_array_finalize() or
 * sync_array_next_idx().
 */
size_t sync_array_available(sync_array_t *sa);

/**
 * Get the index of the next element to become available.
 *
 * This is the same as the number of elements currently in the array.
 *
 * @param sa Array.
 * @return Next element index.
 */
size_t sync_array_next_idx(sync_array_t *sa);

/**
 * Wait for a given element (or any successors) to become available.
 *
 * The waited-for element is guaranteed to be accessible, however its
 * precise location in memory is not guaranteed.  Therefore you must
 * use sync_array_get() to retrieve it.
 *
 * @param sa Array.
 * @param idx Index in the array.
 * @param sec Seconds in timeout, or -1 to wait forever.
 * @param nsec Nanoseconds in timeout.
 * @return 0 for success, <0 for timeout or error.
 */
int sync_array_wait(sync_array_t *sa, size_t idx, int sec, int nsec);

/**
 * Get an element from the array.
 *
 * This function does bounds-checking.  It is not possible to get an
 * element that has been released by all consumers, nor is it possible
 * to get an element that is not yet available.  Best practice is for
 * consumers not to try to get elements that they have released using
 * sync_array_release(), or elements they have not waited for using
 * sync_array_wait().
 *
 * @param sa Array.
 * @param idx Index in the array.
 * @param ent Memory allocated to receive the element.
 * @return 0 for success, <0 for error.
 */
int sync_array_get(sync_array_t *sa, size_t idx, void *out_ent);

/**
 * Append an element to the array.
 *
 * It is forbidden to append elements to a finalized array.
 *
 * @param sa Array.
 * @param ent Element (will be copied).
 * @return 0 for success, <0 for error.
 */
int sync_array_append(sync_array_t *sa, void *ent);

/**
 * Finalize an array.
 *
 * Once an array has been finalized, no more elements can be added to
 * it.  In addition, any consumers waiting for elements after this
 * final index will be resumed, with sync_array_wait() returning a
 * timeout (even if an infinite wait was requested).
 *
 * @param sa Array.
 * @return Final number of elements in the array.
 */
size_t sync_array_finalize(sync_array_t *sa);

/**
 * Forcibly finalize an array.
 *
 * After forcibly finalizing an array, it is no longer possible to add
 * or wait for any elements.  Any consumers waiting will be resumed,
 * with sync_array_wait() returning a timeout (even if an infinite
 * wait was requested), and subsequent calls to sync_array_wait() will
 * return an error.  It is, of course, up to those consumers to exit
 * cleanly.
 *
 * It is recommended that consumers release all elements of the array
 * by calling sync_array_release_all() upon encountering such an error.
 *
 * @param sa Array.
 * @return 0, or <0 on failure (which is not allowed to happen)
 */
int sync_array_force_quit(sync_array_t *sa);

/**
 * Reset an array.
 *
 * This empties the array (i.e. it sets the next index to zero) and
 * unfinalizes it.  It is the responsibility of the caller to notify
 * any consumers of this fact.
 *
 * @param sa Array.
 * @return 0, or <0 on failure (which is not allowed to happen)
 */
int sync_array_reset(sync_array_t *sa);

/**
 * Release elements in an array.
 *
 * sync_array_t can be used as a queue.  In this case, consumers
 * release elements after they have been processed.  When all
 * consumers release all elements before an index, those elements will
 * be freed and their memory returned to the operating system.
 *
 * @param sa Array.
 * @param start_idx First index to release.
 * @param end_idx One past the last index to release.
 * @return Actual last index released.
 */
size_t sync_array_release(sync_array_t *sa,
			  size_t start_idx, size_t end_idx);

/**
 * Release all elements in an array.
 *
 * This fuction is called by consumers to release all elements in the
 * array, such as when they wish to signal completion or respond to a
 * timeout error.
 *
 * @param sa Array.
 * @return Actual last index released.
 */
size_t sync_array_release_all(sync_array_t *sa);

#endif /* __SYNC_ARRAY_H__ */
