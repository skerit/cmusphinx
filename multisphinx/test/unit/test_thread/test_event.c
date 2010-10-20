#include <stdio.h>
#include <sbthread.h>
#include <err.h>
#include <time.h>

int
worker_main(sbthread_t *th)
{
	sbevent_t *cond;
	struct timespec ts;

	cond = sbthread_arg(th);
	printf("Thread %p started\n", th);

	/* Get the first signal. */
	sbevent_wait(cond, -1, -1);
	printf("Thread %p got signal\n", th);

	/* Now wait a while and exit. */
	ts.tv_sec = 0;
	ts.tv_nsec = 500 * 1000 * 1000;
	nanosleep(&ts, NULL);

	return 0;
}

int
main(int argc, char *argv[])
{
	sbthread_t *worker[5];
	sbevent_t *cond;
	int i;

	cond = sbevent_init(TRUE);
	for (i = 0; i < 5; ++i) {
		worker[i] = sbthread_start(NULL, worker_main, cond);
	}

	sbevent_signal(cond);

	for (i = 0; i < 5; ++i) {
		sbthread_free(worker[i]);
	}
	sbevent_free(cond);
	return 0;
}
