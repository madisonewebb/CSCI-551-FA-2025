// https://www.codesansar.com/numerical-methods/gauss-seidel-iteration-using-c-programming.htm
//
// Modified by Sam Siewert 11/4/2021 to add more examples and to provide RHS verify.
// 
// Note that for GSIT, the equations must first be organized into diagonally dominant form,
// where the absolute value of the diagonal coefficient is greater than the sum of the absolute
// values of all other coefficients in that same row.  Otherwise, GSIT will not converge.
//
// Note that this example only works for dimension of 3, but could be generalized for any dimension.
// GSIT is normally used for large sparse matrices that are banded around the diagonal - i.e., most of
// the coefficients are on the diagonal or near it.  However, this program provides a nice simple
// example of how GSIT works.
//
// Check work and compare answers to https://www.symbolab.com/solver/system-of-equations-calculator
// Answers also can be checked with MATLAB.
//
// https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
//
// https://youtu.be/ajJD0Df5CsY
// https://youtu.be/lE3nKRpNOdw
//
#include<stdio.h>
#include<math.h>

// Note that this GSIT demonstration only works for specific dimension problems
#define DIM (5)

/* Arrange systems of linear
   equations to be solved in
   diagonally dominant form
   and form equation for each
   unknown and define here
*/

/* Larger system test
   c1 = (50.0/6.0 + c3*(1.0/6.0))
   c2 = c1
   c3 = 160.0/9.0 + c2*(1.0/9.0))
   c4 = c2*(1.0/11.0) + c3*(8.0/11.0) - c4  + c5*(2.0/11.0)
   c5 = c1*(3.0/4.0) + c2*(1.0/4.0)
*/
#define f1(c1,c2,c3,c4,c5)  (50.0 + c3)/6.0
#define f2(c1,c2,c3,c4,c5)  (c1)
#define f3(c1,c2,c3,c4,c5)  (160.0 + c2)/9.0
#define f4(c1,c2,c3,c4,c5)  (c2 + c3*8.0 +2*c5)/11.0
#define f5(c1,c2,c3,c4,c5)  (3.0*c1 + c2)/4.0

// For verification, enter the coefficients into this array to match equations
double a[DIM][DIM] = {{6.0, 0.0, -1.0, 0.0, 0.0}, {-3.0, 3.0, 0.0, 0.0, 0.0}, {0.0, -1.0, 9.0, 0.0, 0.0}, {0.0, -1.0, -8.0, 11.0, -2.0}, {-3.0, -1.0, 0.0, 0.0, 4.0}};
double b[DIM] = {50.0, 0.0, 160.0, 0.0, 0.0};
double x[DIM]={0.0, 0.0, 0.0, 0.0, 0.0};
int n = DIM;


void verify(double a[DIM][DIM], double *b, double *x, int n);
void vector_print(int nr, double *x);


int main(void)
{
    double x0=0, y0=0, z0=0, a0=0, b0=0, x1, y1, z1, a1, b1, e1, e2, e3, e4, e5, e;
    int count=1;

    printf("Enter tolerable error:\n");
    scanf("%lf", &e);

    printf("\nCount\tx\ty\tz\n");

    // if equations are arranged in diagonally dominate form and are not ill-conditioned, GSIT loop
    // should converge.
    //
    // For non-diaonally dominate inputs, GSIT will likely diverge (growing error).
    //
    do
    {
        /* Calculation */
        x1 = f1(x0,y0,z0, a0, b0);
        y1 = f2(x1,y0,z0, a0, b0);
        z1 = f3(x1,y1,z0, a0, b0);
        a1 = f4(x1,y1,z1, a0, b0);
        b1 = f5(x1,y1,z1, a1, b0);

        /* Error */
        e1 = fabs(x0-x1);
        e2 = fabs(y0-y1);
        e3 = fabs(z0-z1);
        e4 = fabs(a0-a1);
        e5 = fabs(b0-b1);

        printf("%d\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%04f, e1=%16.15f, e2=%16.15f, e3=%16.15f, e4=%16.15f, e5=%16.15f\n",count, x1,y1,z1,a1,b1,e1,e2,e3,e4,e5);
        count++;

        /* Set value for next iteration */
        x0 = x1;
        y0 = y1;
        z0 = z1;
        a0 = a1;
        b0 = b1;

    //} while((e1>e) && (e2>e) && (e3>e));
    } while((e1>e) || (e2>e) || (e3>e) || (e4>e) || (e5>e));

    printf("\nGSIT Solution: x=%0.3f, y=%0.3f z = %0.3f, a = %0.3f, b=%03f\n\n",x1,y1,z1,a1,b1);

    // Additional verification steps to assess error in answer.
    //
    x[0] = x1;
    x[1] = y1;
    x[2] = z1;
    x[3] = a1;
    x[4] = b1;

    verify(a, b, x, n);

    return 0;

}


///////////////////////////////////////////////////////////////////
// Name:     vector_print
//
// Purpose:  This function will print out a one-dimensional
//           vector with supplied dimension.
//
// Usage:    vector_print(nr, x);
//
// Input:    nr     - number of rows (must be >= 1)
//           x      - Vector x[nr] to be printed
//
///////////////////////////////////////////////////////////////////
void vector_print(int nr, double *x)
{
    int row_idx;

    for (row_idx = 0; row_idx < nr; row_idx++)
    {
        printf ("%9.4f  \n", x[row_idx]);
    }

    printf("\n");  // Insert a new line at the end
}


///////////////////////////////////////////////////////////////////
// Multiply coefficient matrix "a" by solution vector s to see if
// it matches the expected RHS we started with.
//
//     a   - Matrix a[n][n]
//     b   - Right hand side vector b[n]
//     x   - Computed solution vector
//     n   - Matrix dimensions
////////////////////////////////////////////////////////////////////
void verify(double a[DIM][DIM], double *b, double *x, int n)
{
    int row_idx, col_jdx;
    double rhs[n];

    // for all rows
    for (row_idx=0; row_idx < n; ++row_idx)
    {
        rhs[row_idx] = 0.0;

        // sum up row's column coefficient x solution vector element
        // as we would do for any matrix * vector operation which yields a vector,
        // which should be the RHS
        for (col_jdx=0; col_jdx < n; ++col_jdx)
        {
            rhs[row_idx] += a[row_idx][col_jdx] * x[col_jdx];
        }
    }

    // Compare original RHS "b" to computed RHS
    printf("Computed RHS is:\n");
    vector_print(n, rhs);

    printf("Original RHS is:\n");
    vector_print(n, b);
}

