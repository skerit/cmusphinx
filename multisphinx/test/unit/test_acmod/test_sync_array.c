#include <stdio.h>
#include <unistd.h>
#include <time.h>

#include <sphinxbase/sbthread.h>
#include <sphinxbase/garray.h>
#include <sphinxbase/err.h>
#include <pocketsphinx.h>

#include "pocketsphinx_internal.h"
#include "test_macros.h"

typedef struct sync_array_s {
	int refcount;
	garray_t *data;
	garray_t *count;
	sbmtx_t *mtx;
	sbevent_t *evt;
	size_t final_next_idx;
} sync_array_t;

sync_array_t *
sync_array_init(size_t n_ent, size_t ent_size)
{
	sync_array_t *sa;

	sa = ckd_calloc(1, sizeof(*sa));
	sa->refcount = 1;
	sa->data = garray_init(n_ent, ent_size);
	sa->count = garray_init(n_ent, 1);
	sa->mtx = sbmtx_init();
	sa->evt = sbevent_init();
	sa->final_next_idx = (size_t)-1;

	return sa;
}

sync_array_t *
sync_array_retain(sync_array_t *sa)
{
	sbmtx_lock(sa->mtx);
	if (sa->refcount == 255) {
		sbmtx_unlock(sa->mtx);
		return NULL;
	}
	++sa->refcount;
	sbmtx_unlock(sa->mtx);

	garray_retain(sa->data);
	garray_retain(sa->count);
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
		sbmtx_unlock(sa->mtx);
		return refcount;
	}
	sbmtx_unlock(sa->mtx);
	printf("Freeing sync array\n");
	garray_free(sa->data);
	garray_free(sa->count);
	sbmtx_free(sa->mtx);
	return 0;
}

void *
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
		if (garray_next_idx(sa->data) > idx)
			return garray_void(sa->data, idx);
		if (idx >= sa->final_next_idx)
			return NULL;
		if (nwait > 0)
			return NULL;
		/* If we wait forever here, there's a race condition
		 * between the tests above and the wait here.  So we
		 * poll the array no matter what, and if we had a
		 * timeout, we make sure to only do it once. */
		if (sbevent_wait(sa->evt, tsec, tnsec) < 0)
			return NULL;
		if (sec != -1)
			++nwait;
	}
	/* Never reached. */
	return NULL;
}

void *
sync_array_append(sync_array_t *sa, void *ent)
{
	int zero = 0;
	void *new_ent;

	sbmtx_lock(sa->mtx);
	/* Not allowed to append to a finalized array. */
	if (garray_next_idx(sa->data) == sa->final_next_idx) {
		sbmtx_unlock(sa->mtx);
		return NULL;
	}
	new_ent = garray_append(sa->data, ent);
	garray_append(sa->count, &zero);
	sbevent_signal(sa->evt);
	sbmtx_unlock(sa->mtx);

	return new_ent;
}

int
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

size_t
sync_array_release(sync_array_t *sa, size_t idx)
{
	size_t i;
	int count;

	sbmtx_lock(sa->mtx);
	if (idx < garray_base(sa->count)
	    || idx > garray_next_idx(sa->count)) {
		sbmtx_unlock(sa->mtx);
		return idx;
	}
	/* Increment count for idx. */
	count = ++garray_ent(sa->count, uint8, idx);

	/* Print stuff for debugging. */
	printf("rc %d counts[%d:%d]:", (int)sa->refcount,
	       (int)garray_base(sa->count),
	       (int)garray_next_idx(sa->count));
	for (i = garray_base(sa->count);
	     i < garray_next_idx(sa->count); ++i)
		printf(" %d", (int)garray_ent(sa->count, uint8, i));
	printf("\n");

	/* Release unreachable elements (this assumes that the
	 * producer will retain a pointer to the array). */
	if (count >= sa->refcount - 1) {
		/* FIXME... Assumes everybody releases everything in
		 * order, which might not happen.... */
		printf("Releasing up to %d\n", (int)idx);
		garray_shift_from(sa->count, idx + 1);
		garray_shift_from(sa->data, idx + 1);
		garray_set_base(sa->count, idx + 1);
		garray_set_base(sa->data, idx + 1);
	}
	sbmtx_unlock(sa->mtx);

	return idx;
}

static int
consumer(sbthread_t *th)
{
	sync_array_t *sa = sbthread_arg(th);
	int i;

	printf("Thread %p sleeping %d secs\n",
	       th, sa->refcount);
	sleep(sa->refcount);
	printf("Thread %p woke up\n", th);

	for (i = 0; i < 20; ++i) {
		int *ent = sync_array_wait(sa, i, -1, -1);
		if (ent != NULL) {
			printf("Thread %p got element %d = %d\n",
			       th, i, *ent);
			sync_array_release(sa, i);
		}
	}
	sync_array_free(sa);
	return 0;
}

int
main(int argc, char *argv[])
{
	sync_array_t *sa;
	int i;

	sa = sync_array_init(0, sizeof(int));
	for (i = 0; i < 10; ++i) {
		struct timespec foo;
		sbthread_start(NULL, consumer,
			       sync_array_retain(sa));
		foo.tv_sec = 0;
		foo.tv_nsec = 50000;
		nanosleep(&foo, NULL);
	}

	/* Now just stream a bunch of numbers in there with some
	 * sleeping in between. */
	for (i = 0; i < 20; ++i) {
		printf("Producer appending %d\n", i);
		sync_array_append(sa, &i);
		sleep(1);
	}
	sync_array_free(sa);
	return 0;
}
