//
// Created by Trinity on 2018-10-23.
//

#include "utils.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <pthread.h>

void fill(double *m, int n) {
  int i;
  for (i = 0; i < n * n; ++i)
    m[i] = i;
//    m[i] = 10;
}

void randfill(double *m, int n, int max, int nn) {
  int i, j;
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      m[i * nn + j] = rand() % max;
}
void memset_zone(double *a, int c, int n, int nn) {
  int i;
  for (i = 0; i < n; ++i)
    memset(a + i * nn, c, sizeof(double) * n);
}

void print_mat(const double *m, int n, const char *band, int nn) {
#ifdef PRINT
  int i, j;
  printf("------------------------- %s -------------------------\n", band);
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j)
      printf("%3.0lf\t", m[i * nn + j]);
    printf("\n");
  }
  printf("----------------------------------------\n");
#endif
}

int is_equal_mat(const double *m1, const double *m2, int n, double *mse) {
  int i;
  *mse = 0;
  for (i = 0; i < n * n; ++i) {
    double err = fabs(m1[i] - m2[i]);
    *mse += err * err;
    if (err > 1e-5)
      return 0;
  }
  *mse /= n * n;
  return 1;
}

enum mode_t {
  PRODUCER,
  CONSUMER
};

struct params_t {
  const double *a;
  const double *b;
  const double *standard;
  int n;
  const char *name;
  void (*fun)(const double *, const double *, double *, int, int, int);
  sem_t *sem;
  enum mode_t mode;
  int ntests;
};

void *tworker(void *p) {
  struct params_t *params = (struct params_t *) p;
  char display[100];
  sprintf(display, "PID[%d] n = %d ", getpid(), params->n);
  strcat(display, params->name);

  switch (params->mode) {
    case PRODUCER: {
      TEST_TIME(
          params->fun(params->a, params->b, (double *) params->standard, params->n, params->n, params->n);,
          display
      );
      if (params->sem != NULL)
        for (int i = 0; i < params->ntests; ++i)
          sem_post(params->sem);
      break;
    }
    case CONSUMER: {
      double *c = (double *) calloc(params->n * params->n, sizeof(double));
      TEST_TIME(
          params->fun(params->a, params->b, c, params->n, params->n, params->n);,
          display
      );
      if (params->sem != NULL)
        sem_wait(params->sem);
      double mse;
      if (!is_equal_mat(params->standard, c, params->n, &mse)) {
        printf("%s ERROR! mse = %lf\n", display, mse);
#ifdef PRINT_ERROR
        print_mat(standard, n, "STANDARD", n);
    print_mat(standard, n, "C", n);
#endif
      }
      free(c);
      break;
    }
  }
  free(params);
  return NULL;
}

pthread_t test_std(const double *a,
                   const double *b,
                   double *standard,
                   int n,
                   void (*fun)(const double *, const double *, double *, int, int, int),
                   const char *name,
                   sem_t *slots,
                   int ntests) {
  struct params_t *pp = (struct params_t *) calloc(1, sizeof(struct params_t));
  pp->a = a;
  pp->b = b;
  pp->standard = standard;
  pp->n = n;
  pp->name = name;
  pp->fun = fun;
  pp->sem = slots;
  pp->mode = PRODUCER;
  pp->ntests = ntests;
  if (slots != NULL) {
    pthread_t t;
    pthread_create(&t, NULL, tworker, pp);
    return t;
  } else {
    tworker(pp);
    return 0;
  }
}

pthread_t test_fun(const double *a,
                   const double *b,
                   const double *standard,
                   int n,
                   void (*fun)(const double *, const double *, double *, int, int, int),
                   const char *name,
                   sem_t *slots) {

  struct params_t *pp = (struct params_t *) calloc(1, sizeof(struct params_t));
  pp->a = a;
  pp->b = b;
  pp->standard = standard;
  pp->n = n;
  pp->name = name;
  pp->fun = fun;
  pp->sem = slots;
  pp->mode = CONSUMER;
  pp->ntests = 0;
  if (slots != NULL) {
    pthread_t t;
    pthread_create(&t, NULL, tworker, pp);
    return t;
  } else {
    tworker(pp);
    return 0;
  }
}

void test_std_mp(const double *a,
                 const double *b,
                 double *standard,
                 int n,
                 void (*fun)(const double *, const double *, double *, int, int, int),
                 const char *name,
                 sem_t *slots,
                 int ntests) {

  int pid;
  if (slots == NULL)
    pid = 0;
  else
    pid = fork();
  if (pid == 0) {
    char display[100];
    sprintf(display, "PID[%d] n = %d ", getpid(), n);
    strcat(display, name);
    TEST_TIME(
        fun(a, b, standard, n, n, n);,
        display
    );
    if (slots != NULL) {
      for (int i = 0; i < ntests; ++i)
        sem_post(slots);
      exit(0);
    } else
      return;
  }
}

void test_fun_mp(const double *a,
                 const double *b,
                 const double *standard,
                 int n,
                 void (*fun)(const double *, const double *, double *, int, int, int),
                 const char *name,
                 sem_t *slots) {

  int pid;
  if (slots == NULL)
    pid = 0;
  else
    pid = fork();
  if (pid == 0) {
    double *c = (double *) calloc(n * n, sizeof(double));
    char display[100];
    sprintf(display, "PID[%d] n = %d ", getpid(), n);
    strcat(display, name);
    TEST_TIME(
        fun(a, b, c, n, n, n);,
        display
    );
    if (slots != NULL)
      sem_wait(slots);
    double mse;
    if (!is_equal_mat(standard, c, n, &mse)) {
      printf("%s ERROR! mse = %lf\n", display, mse);
#ifdef PRINT_ERROR
      print_mat(standard, n, "STANDARD", n);
      print_mat(standard, n, "C", n);
#endif
    }
    free(c);
    if (slots == NULL)
      return;
    else
      exit(0);
  }
}

static inline void get_smat_path(int i, char *path) {
  const char *base = "/matmul_demo.smat";
  char n_str[8];
  sprintf(n_str, "%d", i);
  strcpy(path, base);
  strcat(path, n_str);
}

void *get_shm(const char *location, size_t size, int create) {
  int oflag = O_RDWR;
  if (create) {
    oflag |= O_CREAT;
    shm_unlink(location);
  }
  int fd;
  fd = shm_open(location, oflag, 0700);
  if (fd < 0) {
    printf("shm <%s> failed\n", location);
    perror("why");
    return NULL;
  }

  if (create)
    ftruncate(fd, size);
  void *p = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  close(fd);
  return p;
}

void free_shared(int idx) {
  char path[100];
  get_smat_path(idx, path);
  shm_unlink(path);
}

void *alloc_shared(int idx, size_t n) {
  char path[100];
  get_smat_path(idx, path);
  return get_shm(path, n, 1);
}

