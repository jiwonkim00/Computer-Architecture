#include <stdlib.h>
#include <immintrin.h>
#include "matmul.h"

void transpose(const int* inputB, int* transB, const int K, const int N){
    for(int k=0; k<K; k++){
        for(int j=0; j<N; j++){
            transB[j*K+k] = inputB[k*N+j];
        }
    }
}

//given implementation
/* void matmul(const int M, const int N, const int K, const int *inputA, const int *inputB, int *output) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < K; k++) {
                output[i * N + j] += inputA[i * K + k] * inputB[k * N + j];
            }
        }
    }
} */

//transpose
/* void matmul(const int M, const int N, const int K, const int *inputA, const int *inputB, int *output) {

    int* transB = malloc(N*K*sizeof(int));

    transpose(inputB, transB, K, N);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < K; k++) {
                output[i * N + j] += inputA[i * K + k] * transB[j * K + k];
            }
        }
    }
    free(transB);
}
 */

/* //transpose unrolled
void matmul(const int M, const int N, const int K, const int *inputA, const int *inputB, int *output) {

    int* transB = malloc(N*K*sizeof(int));

    transpose(inputB, transB, K, N);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < K; k=k+8) {
                output[i * N + j] += inputA[i * K + k] * transB[j * K + k];
                output[i * N + j] += inputA[i * K + k+1] * transB[j * K + k+1];
                output[i * N + j] += inputA[i * K + k+2] * transB[j * K + k+2];
                output[i * N + j] += inputA[i * K + k+3] * transB[j * K + k+3];
                output[i * N + j] += inputA[i * K + k+4] * transB[j * K + k+4];
                output[i * N + j] += inputA[i * K + k+5] * transB[j * K + k+5];
                output[i * N + j] += inputA[i * K + k+6] * transB[j * K + k+6];
                output[i * N + j] += inputA[i * K + k+7] * transB[j * K + k+7];
            }
        }
    }
    free(transB);
} */

/* //blocking + loop unrolling
void matmul(const int M, const int N, const int K, const int *inputA, const int *inputB, int *output) {

    int* transB = malloc(N*K*sizeof(int));

    transpose(inputB, transB, K, N);

    int bsize = 32;
    int bsize_k = 256;
    for(int ii = 0; ii < M; ii=ii+bsize){
        for(int jj=0; jj < N; jj=jj+bsize){
            for(int kk = 0; kk < K; kk = kk+bsize_k){
                for(int i=ii; i<ii+bsize; i++){
                    for(int j = jj; j<jj+bsize; j++){
                        int sum = 0;
                        for(int k = kk; k<kk+bsize_k; k=k+8){
                            sum += inputA[i * K + k] * transB[j * K + k];
                            sum += inputA[i * K + k+1] * transB[j * K + k+1];
                            sum += inputA[i * K + k+2] * transB[j * K + k+2];
                            sum += inputA[i * K + k+3] * transB[j * K + k+3];
                            sum += inputA[i * K + k+4] * transB[j * K + k+4];
                            sum += inputA[i * K + k+5] * transB[j * K + k+5];
                            sum += inputA[i * K + k+6] * transB[j * K + k+6];
                            sum += inputA[i * K + k+7] * transB[j * K + k+7];
                        }
                        output[i * N + j] += sum;
                    }
                }
            }
        }
    }
    free(transB);
} */

//AVX_512
void matmul(const int M, const int N, const int K, const int *inputA, const int *inputB, int *output) {

    int* transB = malloc(N*K*sizeof(int));

    __m512i temp_A;
    __m512i temp_B;

    __m512i sum1;
    __m512i sum2;

    __mmask8 mask_all = 0xFF;

    int* temp_sum = malloc(64);

    transpose(inputB, transB, K, N);

    int bsize = 32;
    int bsize_k = 256;
    for(int ii = 0; ii < M; ii=ii+bsize){
        for(int jj=0; jj < N; jj=jj+bsize){
            for(int kk = 0; kk < K; kk = kk+bsize_k){
                for(int i=ii; i<ii+bsize; i++){
                    for(int j = jj; j<jj+bsize; j++){
                        int sum = 0;
                        sum2 = _mm512_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                        for(int k = kk; k<kk+bsize_k; k=k+16){

                            temp_A = _mm512_loadu_si512((void *)(inputA+(i*K+k)));
                            temp_B = _mm512_loadu_si512((void *)(transB+(j*K+k)));

                            sum1 = _mm512_mullo_epi32(temp_A, temp_B);

                            sum2 = _mm512_add_epi32(sum1, sum2);
                        }

                        _mm512_mask_storeu_epi64((__m512i*)temp_sum, mask_all, sum2);

                        sum += temp_sum[0];
                        sum += temp_sum[1];
                        sum += temp_sum[2];
                        sum += temp_sum[3];
                        sum += temp_sum[4];
                        sum += temp_sum[5];
                        sum += temp_sum[6];
                        sum += temp_sum[7];
                        sum += temp_sum[8];
                        sum += temp_sum[9];
                        sum += temp_sum[10];
                        sum += temp_sum[11];
                        sum += temp_sum[12];
                        sum += temp_sum[13];
                        sum += temp_sum[14];
                        sum += temp_sum[15];

                        output[i * N + j] += sum;
                    }
                }
            }
        }
    }

    free(temp_sum);
    free(transB);
}
