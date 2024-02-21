#include "utils.h"
#include "matmul.h"
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <sys/stat.h>        /* For mode_t constants */
#include <fcntl.h>           /* For O_* constants */

int main() {
  srand(time(NULL));
  int min = 32, max = 10000, inc = 32;
  double *a = (double *) calloc(max * max, sizeof(double));
  double *b = (double *) calloc(max * max, sizeof(double));
  double *standard = (double *) calloc(max * max, sizeof(double));
  sem_t *multithread_lock;
//  sem_t sem;
//  sem_init(&sem, 2, 0);
//  multithread_lock = &sem;
  multithread_lock = NULL;  // Disable running test case in parallel
  const int ntests = 12;
  int n;
  pthread_t ts[ntests + 1];
  for (n = min; n <= max; n += inc) {
    printf("Begin test suite, current n = %d\n", n);
    memset(standard, 0, sizeof(double) * n * n);
    randfill(a, n, 10, n);
    randfill(b, n, 10, n);
    pthread_t *tp = ts;
    TEST_TIME(
        *tp++ = test_std(a, b, standard, n, mm1_ijk, "[mm1_ijk]", multithread_lock, ntests);
        *tp++ = test_fun(a, b, standard, n, mm1_ikj, "mm1_ikj", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm1_ikj_v, "mm1_ikj_v", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm1_kji, "mm1_kji", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm2, "mm2", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm2_v, "mm2_v", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm3, "mm3", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm4, "mm4", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm4_cache_block_ikj, "mm4_cache_block_ikj", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm4_cache_block_ijk, "mm4_cache_block_ijk", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm5_vstrassen_ijk, "mm5_vstrassen_ijk", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm5_vstrassen_ikj, "mm5_vstrassen_ikj", multithread_lock);
        *tp++ = test_fun(a, b, standard, n, mm5_vstrassen_ikj_v, "mm5_vstrassen_ikj_v", multithread_lock);
        for (int i = 0; multithread_lock != NULL && i < ntests + 1; ++i)
          pthread_join(ts[i], NULL);,
        "Total test suite"
    );
  }
  free(a);
  free(b);
  free(standard);
  return 0;
}


// Multiprocess test

//int main() {
//  srand(time(NULL));
//  int min = 500, max = 2000, inc = 1;
//  size_t max_size = max*max*sizeof(double);
//  double *a = (double *) calloc(max*max, sizeof(double));
//  double *b = (double *) calloc(max*max, sizeof(double));
//  double *standard = (double *) alloc_shared(3, max_size);
//  sem_t *std_available = sem_open("/matmul_demo.std_available", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
//  sem_init(std_available, 12, 0);
////  sem_t *std_available = NULL; // Single process
//
//  const int ntests = 9;
//  int n;
//  for (n = min; n <= max; n += inc) {
//    printf("Begin test suite, current n = %d\n", n);
//    memset(standard, 0, sizeof(double)*n*n);
//    randfill(a, n, 12, n);
//    randfill(b, n, 12, n);
//    TEST_TIME(
//        test_std_mp(a, b, standard, n, mm1_ijk, "[mm1_ijk]", std_available, ntests);
//        test_fun_mp(a, b, standard, n, mm1_ikj, "mm1_ikj", std_available);
//        test_fun_mp(a, b, standard, n, mm1_kji, "mm1_kji", std_available);
//        test_fun_mp(a, b, standard, n, mm2, "mm2", std_available);
//        test_fun_mp(a, b, standard, n, mm3, "mm3", std_available);
//        test_fun_mp(a, b, standard, n, mm4, "mm4", std_available);
//        test_fun_mp(a, b, standard, n, mm4_cache_block_ikj, "mm4_cache_block_ikj", std_available);
//        test_fun_mp(a, b, standard, n, mm4_cache_block_ijk, "mm4_cache_block_ijk", std_available);
//        test_fun_mp(a, b, standard, n, mm5_vstrassen_ijk, "mm5_vstrassen_ijk", std_available);
//        test_fun_mp(a, b, standard, n, mm5_vstrassen_ikj, "mm5_vstrassen_ikj", std_available);
//        while (wait(0)!=-1);,
//        "Total test suite"
//    );
//  }
//  free(a);
//  free(b);
//  free_shared(3);
//  if (std_available!=NULL)
//    sem_close(std_available);
//  return 0;
//}
