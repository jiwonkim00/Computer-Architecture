//
// Created by Li Yanzhe on 2018-10-23.
//

#ifndef MATMUL_DEMO_UTILS_H
#define MATMUL_DEMO_UTILS_H

#include <semaphore.h>
#include <time.h>
#include <stdio.h>
#include <pthread.h>
//#define PRINT_ERROR
//#define PRINT
#define TEST_TIME(X, NAME)                                                                            \
{                                                                                                     \
  struct timespec stop = {0, 0}, start = {0, 0};                                                      \
  clock_gettime(CLOCK_MONOTONIC, &start);                                                             \
  X                                                                                                   \
  clock_gettime(CLOCK_MONOTONIC, &stop);                                                              \
  printf("%s took %lf s\n", NAME,                                                                     \
  ((double) stop.tv_sec + 1.0e-9 * stop.tv_nsec) - ((double) start.tv_sec + 1.0e-9 * start.tv_nsec)); \
}

void fill(double *m, int n);

void randfill(double *m, int n, int max, int nn);
void memset_zone(double *a, int c, int n, int nn);

//#ifdef PRINT
void print_mat(const double *m, int n, const char *band, int nn);
//#endif

pthread_t test_std(const double *a,
                   const double *b,
                   double *standard,
                   int n,
                   void (*fun)(const double *, const double *, double *, int, int, int),
                   const char *name,
                   sem_t *slots,
                   int ntests);

pthread_t test_fun(const double *a,
                   const double *b,
                   const double *standard,
                   int n,
                   void (*fun)(const double *, const double *, double *, int, int, int),
                   const char *name,
                   sem_t *slots);

void test_std_mp(const double *a,
                 const double *b,
                 double *standard,
                 int n,
                 void (*fun)(const double *, const double *, double *, int, int, int),
                 const char *name,
                 sem_t *slots,
                 int ntests);

void test_fun_mp(const double *a,
                 const double *b,
                 const double *standard,
                 int n,
                 void (*fun)(const double *, const double *, double *, int, int, int),
                 const char *name,
                 sem_t *slots);

void *alloc_shared(int idx, size_t n);

void free_shared(int idx);

#endif //MATMUL_DEMO_UTILS_H
