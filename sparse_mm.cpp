/*
 *   Matrix Market I/O example program
 *
 *   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
 *   and copies it to stdout.  This porgram does nothing useful, but
 *   illustrates common usage of the Matrix Matrix I/O routines.
 *   (See http://math.nist.gov/MatrixMarket for details.)
 *
 *   Usage:  a.out [filename] > output
 *
 *
 *   NOTES:
 *
 *   1) Matrix Market files are always 1-based, i.e. the index of the first
 *      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
 *      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
 *      to files.
 *
 *   2) ANSI C requires one to use the "l" format modifier when reading
 *      double precision floating point numbers in scanf() and
 *      its variants.  For example, use "%lf", "%lg", or "%le"
 *      when reading doubles, otherwise errors will occur.
 */

#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "mmio.c"

void coo_to_csr(int *I, int *J, double *val, int nnz, int n, int **row_ptr, int **col_idx, double **csr_val)
{
    *row_ptr = (int *)malloc((n + 1) * sizeof(int));
    if (*row_ptr == NULL)
    {
        printf("Memory allocation failed for row_ptr\n");
        exit(1);
    }
    *col_idx = (int *)malloc(nnz * sizeof(int));
    if (*col_idx == NULL)
    {
        printf("Memory allocation failed for col_idx\n");
        exit(1);
    }
    *csr_val = (double *)malloc(nnz * sizeof(double));
    if (*csr_val == NULL)
    {
        printf("Memory allocation failed for csr_val\n");
        exit(1);
    }

    // Count number of non-zeros in each row
    for (int i = 0; i < nnz; i++)
    {
        (*row_ptr)[I[i]]++;
    }

    // Calculate row pointers
    (*row_ptr)[n] = nnz;
    for (int i = n - 1; i >= 0; i--)
    {
        (*row_ptr)[i] = (*row_ptr)[i + 1] - (*row_ptr)[i];
    }

    // Fill column indices and values
    for (int i = 0; i < n; i++)
    {
        for (int j = (*row_ptr)[i + 1] - 1; j >= (*row_ptr)[i]; j--)
        {
            (*col_idx)[j] = J[j];
            (*csr_val)[j] = val[j];
        }
    }
}

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i, j, *I, *J;
    double *val;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
        exit(1);

    /* reserve memory for matrices */

    I = (int *)malloc(nz * sizeof(int));
    if (I == NULL)
    {
        printf("Memory allocation failed for I\n");
        exit(1);
    }
    J = (int *)malloc(nz * sizeof(int));
    if (J == NULL)
    {
        printf("Memory allocation failed for J\n");
        exit(1);
    }
    val = (double *)malloc(nz * sizeof(double));
    if (val == NULL)
    {
        printf("Memory allocation failed for val\n");
        exit(1);
    }

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    // Read MM-Matrix
    for (i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--; /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f != stdin)
    {
        fclose(f);
    }

    // Print unchanged MM-Matrix
    // mm_write_banner(stdout, matcode);
    // mm_write_mtx_crd_size(stdout, M, N, nz);
    // for (i = 0; i < nz; i++)
    // {fprintf(stdout, "%d %d %20.19g\n", I[i] + 1, J[i] + 1, val[i]);}
    // printf("\n\n******************************\n\n");

    // Convert to CSR format
    int *row_ptr, *col_idx;
    double *csr_val;
    coo_to_csr(I, J, val, nz, N, &row_ptr, &col_idx, &csr_val);

    // Print the CSR matrix
    // printf("CSR matrix:\n");
    // for (i = 0; i <= N; i++)
    // {
    //    printf("row_ptr[%d] = %d\n", i, row_ptr[i]);
    // }
    // for (i = 0; i < row_ptr[N]; i++)
    // {printf("col_idx[%d] = %d, csr_val[%d] = %.20f\n", i, col_idx[i], i, csr_val[i]);}
    // printf("\n\n******************************\n\n");

    // Initialize the input vector (contains only 1's)
    double *x = (double *)malloc(N * sizeof(double));
    if (x == NULL)
    {
        printf("Memory allocation failed for x\n");
        exit(1);
    }
    for (i = 0; i < N; i++)
    {
        x[i] = 1.0;
    }

    // Allocate memory for the output vector
    double *y = (double *)malloc(N * sizeof(double));
    if (y == NULL)
    {
        printf("Memory allocation failed for y\n");
        exit(1);
    }

    // ***** Perform SpMV *****
    for (int r = 0; r < row_ptr[N]; r++)
    {
        double sum = 0.0;
        for (int i = row_ptr[r]; i < row_ptr[r + 1]; i++)
        {
            if (i + 5 < row_ptr[r + 1])
            {
                // Loads the value of x that's needed in the fifth upcoming iteration
                // prefetch(&x[col_idx[i + 5]], 64);
            }

            sum += val[i] * x[col_idx[i]];
        }
        y[r] += sum;
    }

    // Print the output vector
    printf("Output vector:\n");
    for (i = 0; i < N; i++)
    {
        printf("%.20f\n", y[i]);
    }

    // Free memory
    free(x);
    free(y);

    free(I);
    free(J);
    free(val);

    free(row_ptr);
    free(col_idx);
    free(csr_val);

    return 0;
}