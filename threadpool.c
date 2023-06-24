#include "threadpool.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#if (defined(__linux__) || defined(__APPLE__))
#include <sys/time.h>
#endif


#ifdef _GNU_SOURCE
# undef  _XOPEN_SOURCE
# define _XOPEN_SOURCE 600
# undef  _XOPEN_SOURCE_EXTENDED
# define _XOPEN_SOURCE_EXTENDED 1
# undef  _LARGEFILE64_SOURCE
# define _LARGEFILE64_SOURCE 1
# undef  _BSD_SOURCE
# define _BSD_SOURCE 1
# undef  _SVID_SOURCE
# define _SVID_SOURCE 1
# undef  _ISOC99_SOURCE
# define _ISOC99_SOURCE 1
# undef  _POSIX_SOURCE
# define _POSIX_SOURCE 1
# undef  _POSIX_C_SOURCE
# define _POSIX_C_SOURCE 200112L
# undef  _ATFILE_SOURCE
# define _ATFILE_SOURCE 1
#endif

#if (_WIN32 || _WIN64)
// this block code is written by https://stackoverflow.com/questions/5404277/porting-clock-gettime-to-windows
LARGE_INTEGER getFILETIMEoffset()
{
    SYSTEMTIME s;
    FILETIME f;
    LARGE_INTEGER t;

    s.wYear = 1970;
    s.wMonth = 1;
    s.wDay = 1;
    s.wHour = 0;
    s.wMinute = 0;
    s.wSecond = 0;
    s.wMilliseconds = 0;
    SystemTimeToFileTime(&s, &f);
    t.QuadPart = f.dwHighDateTime;
    t.QuadPart <<= 32;
    t.QuadPart |= f.dwLowDateTime;
    return (t);
}

int clock_gettime(struct timeval* tv)
{
    LARGE_INTEGER           t;
    FILETIME            f;
    double                  microseconds;
    static LARGE_INTEGER    offset;
    static double           frequencyToMicroseconds;
    static int              initialized = 0;
    static BOOL             usePerformanceCounter = 0;

    if (!initialized) {
        LARGE_INTEGER performanceFrequency;
        initialized = 1;
        usePerformanceCounter = QueryPerformanceFrequency(&performanceFrequency);
        if (usePerformanceCounter) {
            QueryPerformanceCounter(&offset);
            frequencyToMicroseconds = (double)performanceFrequency.QuadPart / 1000000.;
        }
        else {
            offset = getFILETIMEoffset();
            frequencyToMicroseconds = 10.;
        }
    }
    if (usePerformanceCounter) QueryPerformanceCounter(&t);
    else {
        GetSystemTimeAsFileTime(&f);
        t.QuadPart = f.dwHighDateTime;
        t.QuadPart <<= 32;
        t.QuadPart |= f.dwLowDateTime;
    }

    t.QuadPart -= offset.QuadPart;
    microseconds = (double)t.QuadPart / frequencyToMicroseconds;
    t.QuadPart = microseconds;
    tv->tv_sec = t.QuadPart / 1000000;
    tv->tv_usec = t.QuadPart % 1000000;
    return (0);
}
#endif

//�������߳�ִ��
void *thread_routine(void *arg)
{
    struct timespec abstime;
    int timeout;
    // printf("thread %d is starting\n", (int)pthread_self());
    threadpool_t *pool = (threadpool_t *)arg;
    while(1)
    {
        timeout = 0;
        //�����̳߳�֮ǰ��Ҫ����
        condition_lock(&pool->ready);
        //����
        pool->idle++;
        //�ȴ������������� ���� �յ��̳߳�����֪ͨ
        while(pool->first == NULL && !pool->quit)
        {
            //�����߳������ȴ�
            // printf("thread %d is waiting\n", (int)pthread_self());
            //��ȡ�ӵ�ǰʱ�䣬�����ϵȴ�ʱ�䣬 ���ý��̵ĳ�ʱ˯��ʱ��
#if (_WIN32 || _WIN64)
            clock_gettime(&abstime);
#else
            clock_gettime(CLOCK_REALTIME, &abstime);
#endif
            abstime.tv_sec += 10;
            int status;
            status = condition_timedwait(&pool->ready, &abstime);  //�ú�������������������̷߳��ʣ���������ʱ������
            if(status == ETIMEDOUT)
            {
#if (_WIN32 || _WIN64)
                printf("thread %d wait timed out\n", pthread_self().x);
#else
                printf("thread %d wait timed out\n", (int)pthread_self());
#endif
                timeout = 1;
                break;
            }
        }

        pool->idle--;
        if(pool->first != NULL)
        {
            //ȡ���ȴ�������ǰ�������Ƴ����񣬲�ִ������
            task_t *t = pool->first;
            pool->first = t->next;
            //��������ִ����Ҫ����ʱ�䣬�Ƚ����������̷߳����̳߳�
            condition_unlock(&pool->ready);
            //ִ������
            t->run(t->arg);
            //ִ���������ͷ��ڴ�
            free(t);
            //���¼���
            condition_lock(&pool->ready);
        }

        //�˳��̳߳�
        if(pool->quit && pool->first == NULL)
        {
            pool->counter--;//��ǰ�������߳���-1
            //���̳߳���û���̣߳�֪ͨ�ȴ��̣߳����̣߳�ȫ�������Ѿ����
            if(pool->counter == 0)
            {
                condition_signal(&pool->ready);
            }
            condition_unlock(&pool->ready);
            break;
        }
        //��ʱ�����������߳�
        if(timeout == 1)
        {
            pool->counter--;//��ǰ�������߳���-1
            condition_unlock(&pool->ready);
            break;
        }

        condition_unlock(&pool->ready);
    }

    // printf("thread %d is exiting\n", (int)pthread_self());
    return NULL;

}


//�̳߳س�ʼ��
void threadpool_init(threadpool_t *pool, int threads)
{

    condition_init(&pool->ready);
    pool->first = NULL;
    pool->last =NULL;
    pool->counter =0;
    pool->idle =0;
    pool->max_threads = threads;
    pool->quit =0;

}


//����һ�������̳߳�
void threadpool_add_task(threadpool_t *pool, void *(*run)(void *arg), void *arg)
{
    //����һ���µ�����
    task_t *newtask = (task_t *)malloc(sizeof(task_t));
    newtask->run = run;
    newtask->arg = arg;
    newtask->next=NULL;//�¼ӵ�������ڶ���β��

    //�̳߳ص�״̬������̹߳�������ǰ��Ҫ����
    condition_lock(&pool->ready);

    if(pool->first == NULL)//��һ���������
    {
        pool->first = newtask;
    }
    else
    {
        pool->last->next = newtask;
    }
    pool->last = newtask;  //����βָ���¼�����߳�

    //�̳߳������߳̿��У�����
    if(pool->idle > 0)
    {
        condition_signal(&pool->ready);
    }
    //��ǰ�̳߳����̸߳���û�дﵽ�趨�����ֵ������һ���µ��߳�
    else if(pool->counter < pool->max_threads)
    {
        pthread_t tid;
        pthread_create(&tid, NULL, thread_routine, pool);
        pool->counter++;
    }
    //����������
    condition_unlock(&pool->ready);
}

//�̳߳�����
void threadpool_destroy(threadpool_t *pool)
{
    //����Ѿ��������٣�ֱ�ӷ���
    if(pool->quit)
    {
    return;
    }
    //����
    condition_lock(&pool->ready);
    //�������ٱ��Ϊ1
    pool->quit = 1;
    //�̳߳����̸߳�������0
    if(pool->counter > 0)
    {
        //���ڵȴ����̣߳������źŻ���
        if(pool->idle > 0)
        {
            condition_broadcast(&pool->ready);
        }
        //����ִ��������̣߳��ȴ����ǽ�������
        while(pool->counter)
        {
            condition_wait(&pool->ready);
        }
    }
    condition_unlock(&pool->ready);
    condition_destroy(&pool->ready);
}