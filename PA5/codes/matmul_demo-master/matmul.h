//
// Created by Trinity on 2018-10-23.
//

#ifndef MATMUL_DEMO_MATMUL_H
#define MATMUL_DEMO_MATMUL_H

void mm1_ijk(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm2_v(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm1_ikj(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm1_kji(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm1_ikj_v(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm2(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm3(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm4(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm4_cache_block_ijk(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm4_cache_block_ikj(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm5_vstrassen_ijk(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm5_vstrassen_ikj(const double *a, const double *b, double *c, int n1, int n2, int n3);

void mm5_vstrassen_ikj_v(const double *a, const double *b, double *c, int n1, int n2, int n3);

//void mmm7(const double *a, const double *b, double *standard, int n1, int n2, int n3);

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
);

#endif //MATMUL_DEMO_MATMUL_H
