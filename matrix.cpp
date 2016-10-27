#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#define VARIABLE_ARRAY		//! 支持可变长数组

#define FABS(x)		mathop(fabs)(x)

real_t matrix_get(matrix_t *m, unsigned int row, unsigned int col)
{
    return MATRIX(m, row, col);
}

void matrix_set(matrix_t *m, unsigned int row, unsigned int col, real_t x)
{
    MATRIX(m, row, col) = x;
}

void matrix_map(matrix_t *s, matrix_t *x, real_t (*f)(real_t))
{
    unsigned int i, j;

    LM_ASSERT(a->row == b->row);
    LM_ASSERT(a->col == b->col);

    matrix_reshape(s, x->row, x->col);

    for(i=0; i<x->row; i++)
        for(j=0; j<x->col; j++)
            MATRIX(s, i, j) = f(MATRIX(x, i, j));
}

void matrix_swap_row(matrix_t *m, unsigned int i, unsigned int j)
{
    unsigned int k;
    real_t tmp;

    LM_ASSERT(i < m->row);
    LM_ASSERT(j < m->row);

    for(k=0; k<m->col; k++)
    {
        tmp = MATRIX(m, i, k);
        MATRIX(m, i, k) = MATRIX(m, j, k);
        MATRIX(m, j, k) = tmp;
    }
}

void matrix_swap_col(matrix_t *m, unsigned int i, unsigned int j)
{
    unsigned int k;
    real_t tmp;

    LM_ASSERT(i < m->col);
    LM_ASSERT(j < m->col);

    for(k=0; k<m->row; k++)
    {
        tmp = MATRIX(m, k, i);
        MATRIX(m, k, i) = MATRIX(m, k, j);
        MATRIX(m, k, j) = tmp;
    }
}

void matrix_copy(matrix_t *to, matrix_t *from)
{
    unsigned int i, j;

    matrix_reshape(to, from->row, from->col);

    for(i=0; i<from->row; i++)
        for(j=0; j<from->col; j++)
            MATRIX(to, i, j) = MATRIX(from, i, j);
}

void matrix_transpose(matrix_t *to, matrix_t *from)
{
    unsigned int i, j;

    matrix_reshape(to, from->col, from->row);

    for(i=0; i<from->row; i++)
        for(j=0; j<from->col; j++)
            MATRIX(to, j, i) = MATRIX(from, i, j);
}

void matrix_eye(matrix_t *m, unsigned int order)
{
    int i, j;

    matrix_reshape(m, order, order);

    for(i=0; i<order; i++)
        for(j=0; j<order; j++)
            MATRIX(m, i, j) = 0;

    for(i=0; i<order; i++)
        MATRIX(m, i, i) = 1;
}

void matrix_add(matrix_t *sum, matrix_t *a, matrix_t *b)
{
    unsigned int i, j;

    LM_ASSERT(a->row == b->row);
    LM_ASSERT(a->col == b->col);

    matrix_reshape(sum, a->row, a->col);

    for(i=0; i<a->row; i++)
        for(j=0; j<a->col; j++)
            MATRIX(sum, i, j) = MATRIX(a, i, j) + MATRIX(b, i, j);
}

void matrix_sub(matrix_t *sub, matrix_t *a, matrix_t *b)
{
    unsigned int i, j;

    LM_ASSERT(a->row == b->row);
    LM_ASSERT(a->col == b->col);

    matrix_reshape(sub, a->row, a->col);

    for(i=0; i<a->row; i++)
        for(j=0; j<a->col; j++)
            MATRIX(sub, i, j) = MATRIX(a, i, j) - MATRIX(b, i, j);
}

void matrix_mul_(matrix_t *mul, matrix_t *x, real_t k)
{
    unsigned int i, j;

    LM_ASSERT(a->row == b->row);
    LM_ASSERT(a->col == b->col);

    matrix_reshape(mul, x->row, x->col);

    for(i=0; i<x->row; i++)
        for(j=0; j<x->col; j++)
            MATRIX(mul, i, j) = k * MATRIX(x, i, j);
}

void matrix_mul(matrix_t *mul, matrix_t *a, matrix_t *b)
{
    unsigned int i, j;

    LM_ASSERT(a->row == b->row);
    LM_ASSERT(a->col == b->col);

    matrix_reshape(mul, a->row, a->col);

    for(i=0; i<a->row; i++)
        for(j=0; j<a->col; j++)
            MATRIX(mul, i, j) = MATRIX(a, i, j) * MATRIX(b, i, j);
}

void matrix_div(matrix_t *div, matrix_t *a, matrix_t *b)
{
    unsigned int i, j;

    LM_ASSERT(a->row == b->row);
    LM_ASSERT(a->col == b->col);

    matrix_reshape(div, a->row, a->col);

    for(i=0; i<a->row; i++)
        for(j=0; j<a->col; j++)
            MATRIX(div, i, j) = MATRIX(a, i, j) / MATRIX(b, i, j);
}

void matrix_pow(matrix_t *p, matrix_t *x, real_t y)
{
    unsigned int i, j;

    LM_ASSERT(a->row == b->row);
    LM_ASSERT(a->col == b->col);

    matrix_reshape(p, x->row, x->col);

    for(i=0; i<x->row; i++)
        for(j=0; j<x->col; j++)
            MATRIX(p, i, j) = mathop(pow)(MATRIX(x, i, j), y);
}

void matrix_sqrt(matrix_t *s, matrix_t *x)
{
    matrix_map(s, x, mathop(sqrt));
}

/**
 * @brief 自然对数
 */
void matrix_ln(matrix_t *s, matrix_t *x)
{
    matrix_map(s, x, mathop(log));
}

/**
 * @brief 以10为底对数
 */
void matrix_log(matrix_t *s, matrix_t *x)
{
    matrix_map(s, x, mathop(log10));
}

/**
 * @brief e的x次方
 */
void matrix_exp(matrix_t *s, matrix_t *x)
{
    matrix_map(s, x, mathop(exp));
}

void matrix_sin(matrix_t *s, matrix_t *x)
{
    matrix_map(s, x, mathop(sin));
}

void matrix_cos(matrix_t *s, matrix_t *x)
{
    matrix_map(s, x, mathop(cos));
}

void matrix_tan(matrix_t *s, matrix_t *x)
{
    matrix_map(s, x, mathop(tan));
}

void matrix_matmul(matrix_t *mul, matrix_t *a, matrix_t *b)
{
    unsigned int i, j, k;

    LM_ASSERT(a->col == b->row);

    matrix_reshape(mul, a->row, b->col);

    for(i=0; i<mul->row; i++)
    {
        for(j=0; j<mul->col; j++)
        {
            MATRIX(mul, i, j) = 0;

            for(k=0; k<a->col; k++)
            {
                MATRIX(mul, i, j) += MATRIX(a, i, k) * MATRIX(b, k, j);
            }
        }
    }
}

int matrix_inv(matrix_t *inv, matrix_t *a)
{
    int i, j, k;
    int ret = 0;

    //! 必须是方阵
    LM_ASSERT(a->row == a->col);

#ifdef VARIABLE_ARRAY
    unsigned int is[a->row];
    unsigned int js[a->col];
#else
    unsigned int *is = malloc(a->row * sizeof(unsigned int));
    unsigned int *js = malloc(a->col * sizeof(unsigned int));
#endif
    real_t max;

    matrix_reshape(inv, a->row, a->col);

    matrix_copy(inv, a);

    for(k=0; k<inv->row; k++)
    {
        //step 1, 全选主元
        max = 0;
        is[k] = k;
        js[k] = k;

        for(i=k; i<inv->row; i++)
        {
            for(j=k; j<inv->col; j++)
            {
                if(max < FABS(MATRIX(inv, i, j)))
                {
                    max = FABS(MATRIX(inv, i, j));
                    is[k] = i;
                    js[k] = j;
                }
            }
        }

        if(max < 0.0001)
        {	//! 无逆矩阵
            ret = -1;
            goto end;
        }

        //交换
        if(is[k] != k)
        {
            matrix_swap_row(inv, k, is[k]);
        }
        if(js[k] != k)
        {
            matrix_swap_col(inv, k, js[k]);
        }

        MATRIX(inv, k, k) = 1 / MATRIX(inv, k, k);

        for(j=0; j<inv->col; j++)
        {
            if(j != k)
                MATRIX(inv, k, j) *= MATRIX(inv, k, k);
        }
        for(i=0; i<inv->row; i++)
        {
            if(i != k)
            {
                for(j=0; j<inv->col; j++)
                {
                    if(j != k)
                        MATRIX(inv, i, j) -= MATRIX(inv, i, k) * MATRIX(inv, k, j);
                }
            }
        }
        for(i=0; i<inv->row; i++)
        {
            if(i != k)
                MATRIX(inv, i, k) *= -MATRIX(inv, k, k);
        }

    }

    //恢复
    //本来 row <-> is[k], column <-> js[k]
    //恢复时：row <-> js[k], column <-> is[k]
    for(k=inv->row-1; k>=0; k--)
    {
        if(js[k] != k)
        {
            matrix_swap_row(inv, k, js[k]);
        }
        if(is[k] != k)
        {
            matrix_swap_col(inv, k, is[k]);
        }
    }

end:
#ifndef VARIABLE_ARRAY
    free(is);
    free(js);
#endif
    return ret;;
}

void matrix_init(matrix_t *m)
{
    LM_ASSERT(m != NULL);

    m->m = NULL;
    m->row = m->col = 0;
}

void matrix_init_with_shape(matrix_t *m, unsigned int row, unsigned int col)
{
    LM_ASSERT(m != NULL);
    LM_ASSERT(row != 0);
    LM_ASSERT(col != 0);

    m->m = (real_t *)malloc(row * col * sizeof(real_t));
    if (m->m != NULL)
    {
        m->row = row;
        m->col = col;
    }
    else
    {
        m->row = m->col = 0;
    }
}

int matrix_reshape(matrix_t *m, unsigned int row, unsigned int col)
{
    LM_ASSERT(m != NULL);
    LM_ASSERT(row != 0);
    LM_ASSERT(col != 0);

    if (row * col == m->row * m->col)
    {
        m->row = row;
        m->col = col;
    }
    else
    {
        if (m->m != NULL)
            free(m->m);

        m->m = (real_t*) malloc(row * col * sizeof(real_t));
        if (m->m != NULL)
        {
            m->row = row;
            m->col = col;
        }
        else
        {
            m->row = m->col = 0;
            return -1;
        }
    }

    return 0;
}

void matrix_release(matrix_t *m)
{
    if(m->m != NULL)
    {
        free(m->m);

        m->m = NULL;
        m->row = m->col = 0;
    }
}

matrix_t *matrix_new(void)
{
    matrix_t *m;

    m = (matrix_t*) malloc(sizeof(matrix_t));
    if (m != NULL)
    {
        matrix_init(m);
    }

    return m;
}

matrix_t *matrix_new_with_shape(unsigned int row, unsigned int col)
{
    matrix_t *m;

    m = (matrix_t*) malloc(sizeof(matrix_t));
    if (m != NULL)
    {
        matrix_init_with_shape(m, row, col);
    }

    return m;
}

void matrix_delete(matrix_t *m)
{
    LM_ASSERT(m != NULL);

    matrix_release(m);
    free(m);
}


real_t matrix_det(matrix *m)
{

}

bool matrix_lu(matrix *L, matrix_t *U, matrix_t *m)
{
    /*This part haven't totally complished! we should use the det of the matrix*/
    int i, j, k;
    real_t s, t;

    int n = m->row;
    cout << "n = " << n << endl;
    matrix *input;
    input = matrix_new();
    matrix_reshape(input, n, n);
    matrix_copy(input, m);

    for(j = 0; j < n; j++)
        matrix_set(input, 0, j, matrix_get(input, 0, j));
    for(i = 1; i < n; i++)
        matrix_set(input, i, 0, matrix_get(input, i, 0) / matrix_get(input, 0, 0) );

    for(k = 1; k < n; k++)
    {
        for(j = k; j < n; j++)
        {
            s = 0;
            for (i = 0; i < k; i++)
                s += matrix_get(input, k, i) * matrix_get(input, i, j); //sum

            matrix_set(input, k, j, matrix_get(input, k, j) - s);
        }
        for(i = k+1; i < n; i++)
        {
            t = 0;
            for(j = 0; j < k; j++)
                t += matrix_get(input, i, j) * matrix_get(input, j, k); //sum
            matrix_set(input, i, k, (matrix_get(input, i, k) - t) / matrix_get(input, k, k));
        }
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(i > j)
            {
                matrix_set(L, i, j, matrix_get(input, i, j));
                matrix_set(U, i, j, 0.);
            }
            else if (i <= j)
            {
                matrix_set(U, i, j, matrix_get(input, i, j));
            }
            if(i == j)
            {
                matrix_set(L, i, j, 1);
            }
        }
    }


}

/*   cholesky decomposition  */
bool matrix_cholesky(matrix_t *L, matrix_t *input)
{
    /*Judge the matrix*/ /*This part haven't totally complished! we should use the det of the matrix*/
    bool flag = true;
    for(int i = 0; i < input->row; i++)
    {
        for(int j = 0; j < input->col; j++)
        {
            if(matrix_get(input, i, j) != matrix_get(input, j, i))
            {
                flag = false;
                cout << "The matrix is not symmetry, Please check out!" << endl;
                break;
            }
            if(i == j)
            {
                if(matrix_get(input, i, j) <= 0)
                {
                    flag = false;
                    cout << "Diagonal elements must be greater than zero(>0)" << endl;
                    break;
                }
            }
        }
    }


    /*cholesky decomposition*/
    if(flag)
    {
        for(int k = 0; k < input->row; k++)
        {
            real_t sum = 0;
            for(int i = 0; i < k; i++)
            {
                sum += matrix_get(L, k, i) * matrix_get(L, k, i);
            }
            sum = matrix_get(input, k, k) - sum;
            matrix_set(L, k, k, sqrt(sum > 0 ? sum: 0));
            for(int i = k+1; i < input->row; i++)
            {
                sum = 0;
                for(int j = 0; j < k; j++)
                {
                    sum += matrix_get(L, i, j) * matrix_get(L, k, j);
                }
                matrix_set(L, i, k, ( (matrix_get(input, i, k) - sum) / matrix_get(L, k, k)));

            }
            for (int j = 0; j < k; j++)
            {
                matrix_set(L, j, k, 0.);
            }
        }
    }

}


/*  Judge the matrix is or not equal to the other one  */
bool matrixIsEqual (matrix_t *input1, matrix_t *input2)
{
    int m = input1->row;
    int n = input1->col;
    int flag = 1;
    int sum = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if ( abs(matrix_get(input1, i, j) - matrix_get(input2, i, j)) <= 0.0001 )
            {
                flag = 1;
            }
            else
            {
                flag = 0;
            }

            sum += flag;

        }
    }
    if (sum == m*n)
        return true;
    else
        return false;

}

bool readFile (const char* filename, const char * checkfile)
{
    /*  Read matrix from files  */
    ifstream inFile;
    ifstream checkFile;
    checkFile.open(checkfile);

    inFile.open(filename);

    clock_t begin,end;
    int sequence;
    int Size = 0;
    float value;
    float checkValue;
    char c;

    double time;

    if(!inFile.is_open())
    {
        cout << "Could not open the file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    if(!checkFile.is_open())
    {
        cout << "Could not open the file " << checkfile << endl;
        exit(EXIT_FAILURE);
    }


    matrix_t *temp, *tempInverse, *tempCheck;

    temp = matrix_new();
    tempInverse = matrix_new();
    tempCheck = matrix_new();

    matrix_reshape(temp, Size, Size);
    matrix_reshape(tempInverse, Size, Size);
    matrix_reshape(tempCheck, Size, Size);

//    matrix_t temp(Size, Size);
//    matrix_t tempInverse(Size, Size);
//    matrix_t tempCheck(Size, Size);

    inFile >> c >> Size >> c >> sequence;

    /*  write the result to files   */
    string outFilename;
    std::stringstream ss;
    ss << Size;
    ss >> outFilename;
    outFilename += "inverse.dat";

    ofstream outFile;
    outFile.open(outFilename.c_str());


    inFile.close();
    inFile.open(filename);

    while(!inFile.eof())
    {
        checkFile >> c >> Size >> c >> sequence;
        inFile >> c >> Size >> c >> sequence;
        //

        outFile << "* " << Size << " " << "* " << sequence << " ";
        cout << "Sequence is " << sequence << " Size is " << Size;

        matrix_reshape(temp, Size, Size);
        matrix_reshape(tempCheck, Size, Size);
//        temp = Matrix(Size, Size);
//        tempCheck = Matrix(Size, Size);

        for(int i = 0; i < temp->row; i++)
        {
            for (int j = 0; j < temp->col; j++)
            {
                inFile >> value;
                matrix_set(temp, i, j, value);

                checkFile >> checkValue;
                matrix_set(tempCheck, i, j, checkValue);
            }
        }

        begin = clock();

        /* solve the inverse */
        matrix_reshape(tempInverse, Size, Size);
        matrix_inv(tempInverse, temp);

        end = clock();
        time = (double) (end - begin) / CLOCKS_PER_SEC;

        cout << "   cost " << time << "s." << endl;


        string result = "True";
        if(matrixIsEqual(tempInverse, tempCheck))
        {
            result = "True";
        }
        else
        {
            result = "False";
        }



        outFile<<setiosflags(ios::fixed)<<setprecision(6);
        for(int i = 0; i < tempInverse->row; i++)
        {
            for (int j = 0; j < tempInverse->col; j++)
            {
                outFile << matrix_get(tempInverse, i, j) << " ";
            }

        }
        outFile << result <<  "\n";

    }

    matrix_delete(temp);
    matrix_delete(tempCheck);
    matrix_delete(tempInverse);


    inFile.close();
    outFile.close();


}
