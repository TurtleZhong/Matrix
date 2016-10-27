#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>

using namespace std;


#ifdef DEBUG
#define LM_ASSERT(x)	do{ if(!(x)) while(1); }while(0)
#else
#define LM_ASSERT(x)
#endif

typedef float real_t;


#if 0		// double
#define mathop(op)		op
#endif

#if 1		// float
#define mathop(op)		op##f
#endif

typedef struct matrix
{
    unsigned int row;
    unsigned int col;
    real_t *m;
}matrix_t;

#define MATRIX(M, r, c) (*((M)->m + (r)*(M)->col + (c)))

void matrix_init(matrix_t *m);
void matrix_init_with_shape(matrix_t *m, unsigned int row, unsigned int col);
int  matrix_reshape(matrix_t *m, unsigned int row, unsigned int col);
void matrix_release(matrix_t *m);
matrix_t *matrix_new(void);
matrix_t *matrix_new_with_shape(unsigned int row, unsigned int col);
void matrix_delete(matrix_t *m);

void matrix_swap_row(matrix_t *m, unsigned int i, unsigned int j);
void matrix_swap_col(matrix_t *m, unsigned int i, unsigned int j);

void matrix_map(matrix_t *s, matrix_t *x, real_t (*f)(real_t));

void matrix_copy(matrix_t *to, matrix_t *from);
void matrix_transpose(matrix_t *to, matrix_t *from);
void matrix_add(matrix_t *sum, matrix_t *a, matrix_t *b);
void matrix_sub(matrix_t *sub, matrix_t *a, matrix_t *b);
void matrix_mul(matrix_t *mul, matrix_t *a, matrix_t *b);
void matrix_matmul(matrix_t *mul, matrix_t *a, matrix_t *b);
int  matrix_inv(matrix_t *inv, matrix_t *a);

real_t matrix_get(matrix_t *m, unsigned int row, unsigned int col);
void matrix_set(matrix_t *m, unsigned int row, unsigned int col, real_t x);
real_t matrix_det(matrix *m);

bool matrix_lu(matrix *L, matrix_t *U, matrix_t *m);
bool matrix_cholesky(matrix_t *L, matrix_t *input);


bool matrixIsEqual (matrix_t *input1, matrix_t *input2);
bool readFile (const char* filename, const char * checkfile);


#endif
