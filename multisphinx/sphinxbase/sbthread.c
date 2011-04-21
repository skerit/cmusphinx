/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights 
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
 * @file sbthread.c
 * @brief Simple portable thread functions
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include <string.h>

#include "sphinxbase/sbthread.h"
#include "sphinxbase/ckd_alloc.h"
#include "sphinxbase/err.h"

/**
 * Semaphore debugging
 */
#define SBTHREAD_SEMDBG
#if defined SBTHREAD_SEMDBG
#define SEMDBG(x) E_INFO x
#else
#define SEMDBG(x)
#endif

/*
 * Platform-specific parts: threads, mutexes, and signals.
 */
#if (defined(_WIN32) || defined(__CYGWIN__)) && !defined(__SYMBIAN32__)
#define _WIN32_WINNT 0x0400
#include <windows.h>

struct sbthread_s {
    cmd_ln_t *config;
    sbthread_main func;
    void *arg;
    HANDLE th;
    DWORD tid;
};

struct sbevent_s {
    HANDLE evt;
};

struct sbmtx_s {
    CRITICAL_SECTION mtx;
};

/* FIXME: Implement sbsem_t for Win32. */

DWORD WINAPI
sbthread_internal_main(LPVOID arg)
{
    sbthread_t *th = (sbthread_t *)arg;
    int rv;

    logfp_index_alloc();
    SEMDBG(("Started thread %p\n", th));
    rv = (*th->func)(th);
    return (DWORD)rv;
}

sbthread_t *
sbthread_start(cmd_ln_t *config, sbthread_main func, void *arg)
{
    sbthread_t *th;

    th = ckd_calloc(1, sizeof(*th));
    th->config = config;
    th->func = func;
    th->arg = arg;
    th->th = CreateThread(NULL, 0, sbthread_internal_main, th, 0, &th->tid);
    if (th->th == NULL) {
        sbthread_free(th);
        return NULL;
    }
    return th;
}

int
sbthread_wait(sbthread_t *th)
{
    DWORD rv, exit;

    /* It has already been joined. */
    if (th->th == NULL)
        return -1;

    rv = WaitForSingleObject(th->th, INFINITE);
    if (rv == WAIT_FAILED) {
        E_ERROR("Failed to join thread: WAIT_FAILED\n");
        return -1;
    }
    GetExitCodeThread(th->th, &exit);
    CloseHandle(th->th);
    th->th = NULL;
    return (int)exit;
}

static DWORD
cond_timed_wait(HANDLE cond, int sec, int nsec)
{
    DWORD rv;
    if (sec == -1) {
        rv = WaitForSingleObject(cond, INFINITE);
    }
    else {
        DWORD ms;

        ms = sec * 1000 + nsec / (1000*1000);
        rv = WaitForSingleObject(cond, ms);
    }
    return rv;
}

/* Silvio Moioli: updated to use Unicode */
sbevent_t *
sbevent_init(int manual_reset)
{
    sbevent_t *evt;

    evt = ckd_calloc(1, sizeof(*evt));
    evt->evt = CreateEventW(NULL, manual_reset, FALSE, NULL);
    if (evt->evt == NULL) {
        ckd_free(evt);
        return NULL;
    }
    return evt;
}

void
sbevent_free(sbevent_t *evt)
{
    CloseHandle(evt->evt);
    ckd_free(evt);
}

int
sbevent_signal(sbevent_t *evt)
{
    return SetEvent(evt->evt) ? 0 : -1;
}

int
sbevent_reset(sbevent_t *evt)
{
    return ResetEvent(evt->evt) ? 0 : -1;
}

int
sbevent_wait(sbevent_t *evt, int sec, int nsec)
{
    DWORD rv;

    rv = cond_timed_wait(evt->evt, sec, nsec);
    return rv;
}

sbmtx_t *
sbmtx_init(void)
{
    sbmtx_t *mtx;

    mtx = ckd_calloc(1, sizeof(*mtx));
    InitializeCriticalSection(&mtx->mtx);
    return mtx;
}

int
sbmtx_trylock(sbmtx_t *mtx)
{
    return TryEnterCriticalSection(&mtx->mtx) ? 0 : -1;
}

int
sbmtx_lock(sbmtx_t *mtx)
{
    EnterCriticalSection(&mtx->mtx);
    return 0;
}

int
sbmtx_unlock(sbmtx_t *mtx)
{
    LeaveCriticalSection(&mtx->mtx);
    return 0;
}

void
sbmtx_free(sbmtx_t *mtx)
{
    DeleteCriticalSection(&mtx->mtx);
    ckd_free(mtx);
}

#else /* POSIX */
#include <pthread.h>
#include <sys/time.h>

struct sbthread_s {
    cmd_ln_t *config;
    sbthread_main func;
    void *arg;
    pthread_t th;
};

struct sbevent_s {
    pthread_mutex_t mtx;
    pthread_cond_t cond;
    /* These are protected by mutexes so they don't need to be atomic
     * types (just sayin' this in case there are any DEC Alpha
     * enthusiasts left in the world...) */
    int16 signalled;
    int16 manual_reset;
};

struct sbmtx_s {
    pthread_mutex_t mtx;
};

struct sbsem_s {
    pthread_mutex_t mtx;
    pthread_cond_t cond;
    char *name;
    int value;
};

static void *
sbthread_internal_main(void *arg)
{
    sbthread_t *th = (sbthread_t *)arg;
    int rv;

    err_set_logfp((FILE *)-1);
    rv = (*th->func)(th);
    return (void *)(long)rv;
}

sbthread_t *
sbthread_start(cmd_ln_t *config, sbthread_main func, void *arg)
{
    sbthread_t *th;
    int rv;

    th = ckd_calloc(1, sizeof(*th));
    th->config = config;
    th->func = func;
    th->arg = arg;
    if ((rv = pthread_create(&th->th, NULL, &sbthread_internal_main, th)) != 0) {
        E_ERROR("Failed to create thread: %d\n", rv);
        sbthread_free(th);
        return NULL;
    }
    return th;
}

int
sbthread_wait(sbthread_t *th)
{
    void *exit;
    int rv;

    /* It has already been joined. */
    if (th->th == (pthread_t)-1)
        return -1;

    rv = pthread_join(th->th, &exit);
    if (rv != 0) {
        E_ERROR("Failed to join thread: %d\n", rv);
        return -1;
    }
    th->th = (pthread_t)-1;
    return (int)(long)exit;
}


static int
cond_timed_wait(pthread_cond_t *cond, pthread_mutex_t *mtx, int sec, int nsec)
{
    int rv;
    if (sec == -1) {
        rv = pthread_cond_wait(cond, mtx);
    }
    else {
        struct timeval now;
        struct timespec end;

        gettimeofday(&now, NULL);
        end.tv_sec = now.tv_sec + sec;
        end.tv_nsec = now.tv_usec * 1000 + nsec;
        if (end.tv_nsec > (1000*1000*1000)) {
            sec += end.tv_nsec / (1000*1000*1000);
            end.tv_nsec = end.tv_nsec % (1000*1000*1000);
        }
        gettimeofday(&now, NULL);
        /* NOTE: returns "non-zero" not "less than zero" */
        rv = pthread_cond_timedwait(cond, mtx, &end);
    }
    return rv;
}

sbevent_t *
sbevent_init(int manual_reset)
{
    sbevent_t *evt;
    int rv;

    evt = ckd_calloc(1, sizeof(*evt));
    if ((rv = pthread_mutex_init(&evt->mtx, NULL)) != 0) {
        E_ERROR("Failed to initialize mutex: %d\n", rv);
        ckd_free(evt);
        return NULL;
    }
    if ((rv = pthread_cond_init(&evt->cond, NULL)) != 0) {
        E_ERROR_SYSTEM("Failed to initialize mutex: %d\n", rv);
        pthread_mutex_destroy(&evt->mtx);
        ckd_free(evt);
        return NULL;
    }
    evt->manual_reset = manual_reset;
    return evt;
}

void
sbevent_free(sbevent_t *evt)
{
    pthread_mutex_destroy(&evt->mtx);
    pthread_cond_destroy(&evt->cond);
    ckd_free(evt);
}

int
sbevent_signal(sbevent_t *evt)
{
    int rv;

    pthread_mutex_lock(&evt->mtx);
    evt->signalled = TRUE;
    rv = pthread_cond_broadcast(&evt->cond);
    pthread_mutex_unlock(&evt->mtx);
    return rv;
}

int
sbevent_reset(sbevent_t *evt)
{
    pthread_mutex_lock(&evt->mtx);
    evt->signalled = FALSE;
    pthread_mutex_unlock(&evt->mtx);

    return 0;
}

int
sbevent_wait(sbevent_t *evt, int sec, int nsec)
{
    int rv = 0;

    /* Lock the mutex before we check its signalled state. */
    pthread_mutex_lock(&evt->mtx);
    /* If it's not signalled, then wait until it is. */
    if (!evt->signalled)
        rv = cond_timed_wait(&evt->cond, &evt->mtx, sec, nsec);
    /* Set its state to unsignalled if we were successful. */
    if (rv == 0 && !evt->manual_reset)
        evt->signalled = FALSE;
    /* And unlock its mutex. */
    pthread_mutex_unlock(&evt->mtx);

    return rv;
}

sbmtx_t *
sbmtx_init(void)
{
    sbmtx_t *mtx;

    mtx = ckd_calloc(1, sizeof(*mtx));
    if (pthread_mutex_init(&mtx->mtx, NULL) != 0) {
        ckd_free(mtx);
        return NULL;
    }
    return mtx;
}

int
sbmtx_trylock(sbmtx_t *mtx)
{
    return pthread_mutex_trylock(&mtx->mtx);
}

int
sbmtx_lock(sbmtx_t *mtx)
{
    return pthread_mutex_lock(&mtx->mtx);
}

int
sbmtx_unlock(sbmtx_t *mtx)
{
    return pthread_mutex_unlock(&mtx->mtx);
}

void
sbmtx_free(sbmtx_t *mtx)
{
    pthread_mutex_destroy(&mtx->mtx);
    ckd_free(mtx);
}

sbsem_t *
sbsem_init(char const *name, int value)
{
    sbsem_t *sem;
    int rv;

    sem = ckd_calloc(1, sizeof(*sem));
    if ((rv = pthread_mutex_init(&sem->mtx, NULL)) != 0) {
        E_ERROR("Failed to initialize mutex: %d\n", rv);
        ckd_free(sem);
        return NULL;
    }
    if ((rv = pthread_cond_init(&sem->cond, NULL)) != 0) {
        E_ERROR_SYSTEM("Failed to initialize mutex: %d\n", rv);
        pthread_mutex_destroy(&sem->mtx);
        ckd_free(sem);
        return NULL;
    }
    sem->name = ckd_salloc(name);
    sem->value = value;
    return sem;
}

void
sbsem_free(sbsem_t *sem)
{
    pthread_mutex_destroy(&sem->mtx);
    pthread_cond_destroy(&sem->cond);
    ckd_free(sem->name);
    ckd_free(sem);
}

int
sbsem_down(sbsem_t *sem, int sec, int nsec)
{
    pthread_mutex_lock(&sem->mtx);
    SEMDBG(("entering sbsem_down(%s),%d\n", sem->name, sem->value));
    while (sem->value <= 0) {
        int rv;
        rv = cond_timed_wait(&sem->cond, &sem->mtx, sec, nsec);
        if (rv != 0) {
            pthread_mutex_unlock(&sem->mtx);
            /* Remember this is a positive error code (why...) */
            return -rv;
        }
    }
    --sem->value;
    SEMDBG(("exiting sbsem_down(%s),%d\n", sem->name, sem->value));
    pthread_mutex_unlock(&sem->mtx);
    return 0;
}

int
sbsem_up(sbsem_t *sem)
{
    int rv = 0;

    pthread_mutex_lock(&sem->mtx);
    if (++sem->value > 0)
        rv = pthread_cond_broadcast(&sem->cond);
    SEMDBG(("sbsem_up(%s),%d\n", sem->name, sem->value));
    pthread_mutex_unlock(&sem->mtx);
    return rv;
}

int
sbsem_set(sbsem_t *sem, int count)
{
    int rv;

    pthread_mutex_lock(&sem->mtx);
    sem->value = count;
    rv = pthread_cond_broadcast(&sem->cond);
    SEMDBG(("sbsem_set(%s),%d\n", sem->name, count));
    pthread_mutex_unlock(&sem->mtx);
    return rv;
}

#endif /* not WIN32 */

cmd_ln_t *
sbthread_config(sbthread_t *th)
{
    return th->config;
}

void *
sbthread_arg(sbthread_t *th)
{
    return th->arg;
}

void
sbthread_free(sbthread_t *th)
{
    sbthread_wait(th);
    ckd_free(th);
}
