//
// Created by Trinity on 2018-10-23.
//

#include "utils.h"
#include "matmul.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <memory.h>

#define BSIZE 32
typedef double vector8 __attribute__ ((vector_size (sizeof(double)*8)));

void mm1_ijk(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k;
  for (i = 0; i < n1; ++i)
    for (j = 0; j < n3; ++j)
      for (k = 0; k < n2; ++k)
        c[i * n1 + j] += a[i * n1 + k] * b[k * n1 + j];
}

void mm1_ikj(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k;
  for (i = 0; i < n1; ++i)
    for (k = 0; k < n2; ++k)
      for (j = 0; j < n3; ++j)
        c[i * n1 + j] += a[i * n1 + k] * b[k * n1 + j];
}

void mm1_ikj_v(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k;
  if (n3 % 8 != 0)
    return;
  for (i = 0; i < n1; ++i) {
    register int ti = i * n1;
    for (k = 0; k < n2; ++k) {
      register int tk = k * n2;
      register double x = a[ti + k];
      vector8 vk = {x, x, x, x, x, x, x, x};
      for (j = 0; j < n3; j += 8) {
        vector8 vj =
            {b[tk + j], b[tk + j + 1], b[tk + j + 2], b[tk + j + 3], b[tk + j + 4], b[tk + j + 5], b[tk + j + 6],
             b[tk + j + 7]};
        vector8 ans = vk * vj;
        c[ti + j] += ans[0];
        c[ti + j + 1] += ans[1];
        c[ti + j + 2] += ans[2];
        c[ti + j + 3] += ans[3];
        c[ti + j + 4] += ans[4];
        c[ti + j + 5] += ans[5];
        c[ti + j + 6] += ans[6];
        c[ti + j + 7] += ans[7];
      }
    }
  }
}

void mm1_kji(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k;
  for (k = 0; k < n2; ++k)
    for (j = 0; j < n3; ++j)
      for (i = 0; i < n1; ++i)
        c[i * n1 + j] += a[i * n1 + k] * b[k * n1 + j];
}

void mm2(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k;
  for (i = 0; i < n1; ++i)
    for (j = 0; j < n3; ++j) {
      register double t = 0;
      for (k = 0; k < n2; ++k) {
        t += a[i * n1 + k] * b[k * n1 + j];
      }
      c[i * n1 + j] = t;
    }
}

void mm2_v(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k;
  if (n2 % 8 != 0)
    return;
  for (i = 0; i < n1; ++i)
    for (j = 0; j < n3; ++j) {
      register double t = 0;
      for (k = 0; k < n2; k += 8) {
        vector8 va = {a[i * n1 + k], a[i * n1 + k + 1], a[i * n1 + k + 2], a[i * n1 + k + 3], a[i * n1 + k + 4],
                      a[i * n1 + k + 5],
                      a[i * n1 + k + 6], a[i * n1 + k + 7]};
        vector8 vb =
            {b[k * n1 + j], b[(k + 1) * n1 + j], b[(k + 2) * n1 + j], b[(k + 3) * n1 + j], b[(k + 4) * n1 + j],
             b[(k + 5) * n1 + j],
             b[(k + 6) * n1 + j], b[(k + 7) * n1 + j]};
        vector8 ans = va * vb;
        t += ans[0] + ans[1] + ans[2] + ans[3] + ans[4] + ans[5] + ans[6] + ans[7];
      }
      c[i * n1 + j] = t;
    }
}

// Requirements:
// - Square only
// - n is even
void mm3(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k;
  if (n1 % 2 != 0 || n2 % 2 != 0 || n3 % 2 != 0)
    return;
  for (i = 0; i < n1; i += 2) {
    for (j = 0; j < n3; j += 2) {
      for (k = 0; k < n2; k += 2) {
        c[i * n1 + j] += a[i * n1 + k] * b[k * n1 + j] + a[i * n1 + k + 1] * b[(k + 1) * n1 + j];
        c[(i + 1) * n1 + j] += a[(i + 1) * n1 + k] * b[k * n1 + j] + a[(i + 1) * n1 + k + 1] * b[(k + 1) * n1 + j];
        c[i * n1 + j + 1] += a[i * n1 + k] * b[k * n1 + j + 1] + a[i * n1 + k + 1] * b[(k + 1) * n1 + j + 1];
        c[(i + 1) * n1 + j + 1] +=
            a[(i + 1) * n1 + k] * b[k * n1 + j + 1] + a[(i + 1) * n1 + k + 1] * b[(k + 1) * n1 + j + 1];
      }
    }
  }
}

// Requirements:
// - Square only
// - n is even
void mm4(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k;
  if (n1 % 2 != 0 || n2 % 2 != 0 || n3 % 2 != 0)
    return;
  for (i = 0; i < n1; i += 2) {
    for (j = 0; j < n3; j += 2) {
      register double c00 = 0, c01 = 0, c10 = 0, c11 = 0;
//      register double c00 = standard[i * n + j], c01 = standard[i * n + j + 1], c10 = standard[(i + 1) * n + j], c11 = standard[(i + 1) * n + j];
      for (k = 0; k < n2; k += 2) {
        register double a00 = a[i * n1 + k], a01 = a[i * n1 + k + 1], a10 = a[(i + 1) * n1 + k],
            a11 = a[(i + 1) * n1 + k + 1];
        register double b00 = b[k * n1 + j], b10 = b[(k + 1) * n1 + j], b01 = b[k * n1 + j + 1],
            b11 = b[(k + 1) * n1 + j + 1];
        c00 += a00 * b00 + a01 * b10;
        c10 += a10 * b00 + a11 * b10;
        c01 += a00 * b01 + a01 * b11;
        c11 += a10 * b01 + a11 * b11;
      }
      c[i * n1 + j] = c00;
      c[i * n1 + j + 1] = c01;
      c[(i + 1) * n1 + j] = c10;
      c[(i + 1) * n1 + j + 1] = c11;
    }
  }
}

// Requirements:
// - Square only
// - n % B == 0
void mm4_cache_block_ijk(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k, i1, j1, k1;
  const int n = n1;
  if (n % BSIZE != 0)
    return;
  for (i = 0; i < n; i += BSIZE)
    for (j = 0; j < n; j += BSIZE)
      for (k = 0; k < n; k += BSIZE)
        for (i1 = i; i1 < i + BSIZE; i1 += 2)
          for (j1 = j; j1 < j + BSIZE; j1 += 2) {
            /* 2 by 2 mini matrix multiplications in registers*/
            register int t = i1 * n + j1;
            register int tt = t + n;
            register double c00 = c[t];
            register double c01 = c[t + 1];
            register double c10 = c[tt];
            register double c11 = c[tt + 1];
            for (k1 = k; k1 < k + BSIZE; k1 += 2) {
              register int ta = i1 * n + k1;
              register int tta = ta + n;
              register int tb = k1 * n + j1;
              register int ttb = tb + n;
              register double b00 = b[tb];
              register double b01 = b[tb + 1];
              register double b10 = b[ttb];
              register double b11 = b[ttb + 1];
              register double a00 = a[ta];
              register double a01 = a[ta + 1];
              register double a10 = a[tta];
              register double a11 = a[tta + 1];
              /*instructions are reordered for better pipeline*/
              c00 += a00 * b00;
              c01 += a00 * b01;
              c10 += a10 * b00;
              c11 += a10 * b01;
              c00 += a01 * b10;
              c01 += a01 * b11;
              c10 += a11 * b10;
              c11 += a11 * b11;
            }
            c[t] = c00;
            c[t + 1] = c01;
            c[tt] = c10;
            c[tt + 1] = c11;
          }
}

// Requirements:
// - Square only
// - n % B == 0
void mm4_cache_block_ikj(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  int i, j, k, i1, j1, k1;
  const int n = n1;
  if (n % BSIZE != 0)
    return;
  for (i = 0; i < n; i += BSIZE)
    for (k = 0; k < n; k += BSIZE)
      for (j = 0; j < n; j += BSIZE)
        for (i1 = i; i1 < i + BSIZE; i1 += 2)
          for (j1 = j; j1 < j + BSIZE; j1 += 2) {
            /* 2 by 2 mini matrix multiplications in registers*/
            register int t = i1 * n + j1;
            register int tt = t + n;
            register double c00 = c[t];
            register double c01 = c[t + 1];
            register double c10 = c[tt];
            register double c11 = c[tt + 1];
            for (k1 = k; k1 < k + BSIZE; k1 += 2) {
              register int ta = i1 * n + k1;
              register int tta = ta + n;
              register int tb = k1 * n + j1;
              register int ttb = tb + n;
              register double b00 = b[tb];
              register double b01 = b[tb + 1];
              register double b10 = b[ttb];
              register double b11 = b[ttb + 1];
              register double a00 = a[ta];
              register double a01 = a[ta + 1];
              register double a10 = a[tta];
              register double a11 = a[tta + 1];
              /*instructions are reordered for better pipeline*/
              c00 += a00 * b00 + a01 * b10;
              c10 += a10 * b00 + a11 * b10;
              c01 += a00 * b01 + a01 * b11;
              c11 += a10 * b01 + a11 * b11;
            }
            c[t] = c00;
            c[t + 1] = c01;
            c[tt] = c10;
            c[tt + 1] = c11;
          }
}

static inline
void scope_mul_ijk(const double *a, const double *b, double *c, int n1, int n2, int n3, int nn) {
  int i, j, k;
  for (i = 0; i < n1; ++i) {
    register int ti = i * nn;
    for (j = 0; j < n3; ++j) {
      register double sum = 0;
      for (k = 0; k < n2; ++k)
        sum += a[ti + k] * b[k * nn + j];
      c[ti + j] += sum;
    }
  }
}
static inline
void scope_mul_ikj(const double *a, const double *b, double *c, int n1, int n2, int n3, int nn) {
  int i, j, k;
  for (i = 0; i < n1; ++i) {
    register int ti = i * nn;
    for (k = 0; k < n2; ++k) {
      register int tk = k * nn;
      for (j = 0; j < n3; ++j)
        c[ti + j] += a[ti + k] * b[tk + j];
    }
  }
}

static inline
void scope_mul_ikj_v(const double *a, const double *b, double *c, int n1, int n2, int n3, int nn) {
  int i, j, k;
  for (i = 0; i < n1; ++i) {
    register int ti = i * nn;
    for (k = 0; k < n2; k += 1) {
      register int tk = k * nn;
      register double x = a[ti + k];
      vector8 vk = {x, x, x, x, x, x, x, x};
      for (j = 0; j < n3; j += 8) {
        if (n3 - j < 8) {
          for (int kj = j; kj < n3; ++kj) {
            c[ti + kj] += x * b[tk + kj];
          }
        } else {
          vector8 vj =
              {b[tk + j], b[tk + j + 1], b[tk + j + 2], b[tk + j + 3], b[tk + j + 4], b[tk + j + 5], b[tk + j + 6],
               b[tk + j + 7]};
          vector8 ans = vk * vj;
          c[ti + j] += ans[0];
          c[ti + j + 1] += ans[1];
          c[ti + j + 2] += ans[2];
          c[ti + j + 3] += ans[3];
          c[ti + j + 4] += ans[4];
          c[ti + j + 5] += ans[5];
          c[ti + j + 6] += ans[6];
          c[ti + j + 7] += ans[7];
        }
      }
    }
  }
}

static inline
void madd_v(const double *a, const double *b, double *c, int m, int n, int nn) {
  int i, j;
  for (i = 0; i < m; ++i) {
    const int head = i * nn;
    for (j = 0; j < n; j += 8) {
      if (n - j < 8) {
        for (int k = j; k < n; ++k) {
          c[head + k] = a[head + k] + b[head + k];
        }
      } else {
        vector8 va = {a[head + j], a[head + j + 1], a[head + j + 2], a[head + j + 3], a[head + j + 4], a[head + j + 5],
                      a[head + j + 6], a[head + j + 7]};
        vector8 vb = {b[head + j], b[head + j + 1], b[head + j + 2], b[head + j + 3], b[head + j + 4], b[head + j + 5],
                      b[head + j + 6], b[head + j + 7]};
        vector8 ans = va + vb;
        c[head + j] = ans[0];
        c[head + j + 1] = ans[1];
        c[head + j + 2] = ans[2];
        c[head + j + 3] = ans[3];
        c[head + j + 4] = ans[4];
        c[head + j + 5] = ans[5];
        c[head + j + 6] = ans[6];
        c[head + j + 7] = ans[7];
      }
    }
  }
}

static inline
void madd(const double *a, const double *b, double *c, int m, int n, int nn) {
  int i, j;
  for (i = 0; i < m; ++i) {
    const int head = i * nn;
    for (j = 0; j < n; ++j)
      c[head + j] = a[head + j] + b[head + j];
  }
}

static inline
void msub_v(const double *a, const double *b, double *c, int m, int n, int nn) {
  int i, j;
  for (i = 0; i < m; ++i) {
    const int head = i * nn;
    for (j = 0; j < n; j += 8) {
      if (n - j < 8) {
        for (int k = j; k < n; ++k) {
          c[head + k] = a[head + k] - b[head + k];
        }
      } else {
        vector8 va = {a[head + j], a[head + j + 1], a[head + j + 2], a[head + j + 3], a[head + j + 4], a[head + j + 5],
                      a[head + j + 6], a[head + j + 7]};
        vector8 vb = {b[head + j], b[head + j + 1], b[head + j + 2], b[head + j + 3], b[head + j + 4], b[head + j + 5],
                      b[head + j + 6], b[head + j + 7]};
        vector8 ans = va - vb;
        c[head + j] = ans[0];
        c[head + j + 1] = ans[1];
        c[head + j + 2] = ans[2];
        c[head + j + 3] = ans[3];
        c[head + j + 4] = ans[4];
        c[head + j + 5] = ans[5];
        c[head + j + 6] = ans[6];
        c[head + j + 7] = ans[7];
      }
    }
  }
}

static inline
void msub(const double *a, const double *b, double *c, int m, int n, int nn) {
  int i, j;
  for (i = 0; i < m; ++i) {
    const int head = i * nn;
    for (j = 0; j < n; ++j)
      c[head + j] = a[head + j] - b[head + j];
  }
}
/*
 [Illustration]

        Draft Zone              NxN Matrix A              NxN Matrix B              NxN Matrix C
      nmx := next mx
  ---------------------     ---------------------     ---------------------     ---------------------
  |         |         |     |         |         |     |         |         |     |         |         |
  |    mx   |   nmx   |     |   a00   |   a01   |     |   b00   |   b01   |     |   c00   |   c01   |
  |         |         |     |         |         |     |         |         |     |         |         |
  ---------------------     ---------------------  X  ---------------------  =  ---------------------
  |         |         |     |         |         |     |         |         |     |         |         |
  |  azone  |  bzone  |     |   a10   |   a11   |     |   b10   |   b11   |     |   c10   |   c11   |
  |         |         |     |         |         |     |         |         |     |         |         |
  ---------------------     ---------------------     ---------------------     ---------------------

  m1 = (a00+a11)*(b00+b11)                            c00 = m1 + m4 - m5 + m7
  m2 = (a10+a11)*b00                                  c01 = m3 + m5
  m3 = a00*(b01-b11)                                  c10 = m2 + m4
  m4 = a11*(b10-b00)                                  c11 = m1 + m3 - m2 + m6
  m5 = (a00+a01)*b11
  m6 = (a10-a00)*(b00+b01)
  m7 = (a01-a11)*(b10+b11)

[Params]

 Inputs:
  a00 := pointer of the element at matrix A's UP-LEFT corner
  b00 := pointer of the element at matrix B's UP-LEFT corner
  n   := the ORDER of current square window of operation
  nn  := the ORDER of original matrix before first recursion, this value shouldn't be changed during recursion
  mx  := the allocated draft zone provided by caller, should be initialized to 0 before first recursion
         length = nn * nn, size = sizeof(double) * nn * nn
 Outputs:
  c00 := pointer of the element at result matrix C's UP-LEFT corner
*/

void v_strassen(
    const double *a00,
    const double *b00,
    double *c00,
    int n,
    int nn,
    double *mx,
    void (*basic_mul)(const double *, const double *, double *, int, int, int, int),
    void (*madd)(const double *, const double *, double *, int, int, int),
    void (*msub)(const double *, const double *, double *, int, int, int)
) {
  if (n <= BSIZE) {
    basic_mul(a00, b00, c00, n, n, n, nn);
  } else if (n % 2 == 0) {
    const int half = n / 2;

    const double *a01 = a00 + half;
    const double *a10 = a00 + half * nn;
    const double *a11 = a10 + half;

    const double *b01 = b00 + half;
    const double *b10 = b00 + half * nn;
    const double *b11 = b10 + half;

    double *c01 = c00 + half;
    double *c10 = c00 + half * nn;
    double *c11 = c10 + half;
    double *nmx = mx + half;
    double *azone = mx + half * nn;
    double *bzone = azone + half;

    memset_zone(mx, 0, half, nn);
    madd(a00, a11, azone, half, half, nn);                              //azone = a00+a11
    madd(b00, b11, bzone, half, half, nn);                              //bzone = b00+b11
    v_strassen(azone, bzone, mx, half, nn, nmx, basic_mul, madd, msub); //m1 = (a00+a11)*(b00+b11)
    madd(c00, mx, c00, half, half, nn);                                 //c00 += m1
    madd(c11, mx, c11, half, half, nn);                                 //c11 += m1

    memset_zone(mx, 0, half, nn);
    madd(a10, a11, azone, half, half, nn);                              //azone = a10+a11
    v_strassen(azone, b00, mx, half, nn, nmx, basic_mul, madd, msub);   //m2 = (a10+a11)*b00
    madd(c10, mx, c10, half, half, nn);                                 //c10 += m2
    msub(c11, mx, c11, half, half, nn);                                 //c11 -= m2

    memset_zone(mx, 0, half, nn);
    msub(b01, b11, bzone, half, half, nn);                              //bzone = b01-b11
    v_strassen(a00, bzone, mx, half, nn, nmx, basic_mul, madd, msub);   //m3 = a00*(b01-b11)
    madd(c01, mx, c01, half, half, nn);                                 //c01 += m3
    madd(c11, mx, c11, half, half, nn);                                 //c11 += m3

    memset_zone(mx, 0, half, nn);
    msub(b10, b00, bzone, half, half, nn);                              //bzone = b10-b00
    v_strassen(a11, bzone, mx, half, nn, nmx, basic_mul, madd, msub);   //m4 = a11*(b10-b00)
    madd(c00, mx, c00, half, half, nn);                                 //c00 += m4
    madd(c10, mx, c10, half, half, nn);                                 //c10 += m4

    memset_zone(mx, 0, half, nn);
    madd(a00, a01, azone, half, half, nn);                              //azone = a00+a01
    v_strassen(azone, b11, mx, half, nn, nmx, basic_mul, madd, msub);   //m5 = (a00+a01)*b11
    msub(c00, mx, c00, half, half, nn);                                 //c00 -= m5
    madd(c01, mx, c01, half, half, nn);                                 //c01 += m5

    msub(a10, a00, azone, half, half, nn);                              //azone = a10-a00
    madd(b00, b01, bzone, half, half, nn);                              //bzone = b00+b01
    v_strassen(azone, bzone, c11, half, nn, nmx, basic_mul, madd, msub);//c11 += (a10-a00)*(b00+b01)

    msub(a01, a11, azone, half, half, nn);                              //azone = a01-a11
    madd(b10, b11, bzone, half, half, nn);                              //bzone = b10+b11
    v_strassen(azone, bzone, c00, half, nn, nmx, basic_mul, madd, msub);//c00 += (a01-a11)*(b10+b11)
  } else {
    const int part = n / 2 + 1;
    const int rest = n - part;
    const double *a01 = a00 + part;
    const double *a10 = a00 + part * nn;
    const double *a11 = a10 + part;

    const double *b01 = b00 + part;
    const double *b10 = b00 + part * nn;
    const double *b11 = b10 + part;

    double *c01 = c00 + part;
    double *c10 = c00 + part * nn;
    double *c11 = c10 + part;
//    double *mx_ur = mx + part;
//    double *mx_bl = mx + part * nn;
//    double *mx_br = mx_bl + part;

    //c00 = a00 * b00 + a01 * b10;
    //c01 = a00 * b01 + a01 * b11;
    //c10 = a10 * b00 + a11 * b10;
    //c11 = a10 * b01 + a11 * b11;
    v_strassen(a00, b00, c00, part, nn, mx, basic_mul, madd, msub);   //c00 += a00 * b00
    basic_mul(a01, b10, c00, part, rest, part, nn);                   //c00 += a01 * b10

    basic_mul(a00, b01, c01, part, part, rest, nn);                   //c01 += a00 * b01
    basic_mul(a01, b11, c01, part, rest, rest, nn);                   //c01 += a01 * b11

    basic_mul(a10, b00, c10, rest, part, part, nn);                   //c10 += a00 * b01
    basic_mul(a11, b10, c10, rest, rest, part, nn);                   //c10 += a11 * b10

    basic_mul(a10, b01, c11, rest, part, rest, nn);                   //c11 += a10 * b01
    v_strassen(a11, b11, c11, rest, nn, mx, basic_mul, madd, msub);   //c11 += a11 * b11
  }
}

void mm5_vstrassen_ijk(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  double *mx = (double *) calloc((size_t) n1 * n1, sizeof(double));
  v_strassen(a, b, c, n1, n1, mx, scope_mul_ijk, madd, msub);
  free(mx);
}

void mm5_vstrassen_ikj(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  double *mx = (double *) calloc((size_t) n1 * n1, sizeof(double));
  v_strassen(a, b, c, n1, n1, mx, scope_mul_ikj, madd, msub);
  free(mx);
}
void mm5_vstrassen_ikj_v(const double *a, const double *b, double *c, int n1, int n2, int n3) {
  double *mx = (double *) calloc((size_t) n1 * n1, sizeof(double));
  v_strassen(a, b, c, n1, n1, mx, scope_mul_ikj_v, madd_v, msub_v);
  free(mx);
}

//void mmm7(const double *a, const double *b, double *c, int n) {
//  int i, j, k;
//  for (i = 0; i < n; i += 2) {
//    for (j = 0; j < n; j += 2) {
//      vector8 r = {0};
//      for (k = 0; k < n; k += 2) {
//        register double a1 = a[i*n + k], a2 = a[i*n + k + 1], a3 = a[(i + 1)*n + k], a4 = a[(i + 1)*n + k + 1];
//        register double b1 = b[k*n + j], b2 = b[(k + 1)*n + j], b3 = b[k*n + j + 1], b4 = b[(k + 1)*n + j + 1];
//        vector8 v1 = {a1, a3, a1, a3};
//        vector8 v2 = {b1, b1, b3, b3};
//        vector8 v3 = {a2, a4, a2, a4};
//        vector8 v4 = {b2, b2, b4, b4};
//        r = v2*v2 + v3*v4 + r;
//      }
//      c[i*n + j] = r[0];
//      c[i*n + j + 1] = r[2];
//      c[(i + 1)*n + j] = r[1];
//      c[(i + 1)*n + j + 1] = r[3];
//    }
//  }
//}
