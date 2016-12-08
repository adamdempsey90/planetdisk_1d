#include "pdisk.h"

extern void zgemm_(char *, char *, int *, int *,int *, double complex *, double complex *, int *, double complex *, int *, double complex *, double complex *, int *);
extern void zgemv_(char *, int *, int *, double complex *, double complex *, int *, double complex *,int *,double complex *,double complex *,int*);
extern void zgesv_(int *, int *, double complex *, int *, int *, double complex *, int *, int *);

void matvec(double *ld, double *md, double *ud, double *y, double *res, double a, double b,int n) {
    /* res = a*res + b*mat.y, 
     * where mat is a tridiagonal matrix with 
     * md = main diagonal
     * ld = lower diagonal
     * ud = upper diagonal
     */
    int i;

    res[0] = a*res[0] + b*(md[0]*y[0] + ud[0]*y[1]);



#ifdef _OPENMP
#pragma omp parallel for private(i) 
#endif
    for(i=1;i<n-1;i++) {
        res[i] = a*res[i] + b*( ld[i-1]*y[i-1] + md[i]*y[i] + ud[i]*y[i+1]);

    }
    res[n-1] = a*res[n-1] + b*( ld[n-2]*y[n-2] + md[n-1] * y[n-1]);

    return;
}

void matvec_full(double *ld, double *md, double *ud,double *u, double *w, double *y, double *res, double a, double b,int n,int nw) {
    /* res = a*res + b*(mat + u1 w1^T + u2 w2^T + ...).y, 
     * where mat is a tridiagonal matrix with 
     * md = main diagonal
     * ld = lower diagonal
     * ud = upper diagonal
     * u and w are vectors which add to mat via an outer product
     */
    int i,j,k;
    double temp;


    res[0] = a*res[0] + b* (md[0]*y[0] + ud[0]*y[1]);

    temp=0;
#ifdef _OPENMP
#pragma omp parallel for private(j,k) collapse(2) reduction(+:temp)  
#endif
    for(j=0;j<nw;j++) {
        for(k=0;k<n;k++) {
           temp += b*u[0+n*j]*w[k + n*j]*y[k]; 
        }
    }
    res[0] += temp;

#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)  
#endif
    for(i=1;i<n-1;i++) {
        res[i] = a*res[i] + b* (ld[i-1]*y[i-1] + md[i]*y[i] + ud[i]*y[i+1]);
        for(j=0;j<nw;j++) {
            for(k=0;k<n;k++) {
                res[i] += b*u[i+n*j]*w[k + n*j]*y[k]; 
            }
        }
    }
    res[n-1] = a*res[n-1] + b*( ld[n-2]*y[n-2] + md[n-1] * y[n-1]);

    temp=0;
#ifdef _OPENMP
#pragma omp parallel for private(j,k) collapse(2) reduction(+:temp)  
#endif
    for(j=0;j<nw;j++) {
        for(k=0;k<n;k++) {
            temp += b*u[n-1+n*j]*w[k + n*j]*y[k]; 
        }
    }
    res[n-1] += temp;
    return;
}
double dotprod(double *v1, double *v2, int n) {
    int i;

    double res=0;

#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(+:res)
#endif
    for(i=0;i<n;i++) {
        res += v1[i]*v2[i];
    }
    return res;


}
void outer(double *v1, double *v2, double *res, int n) {
    int i,j;

#ifdef _OPENMP
#pragma omp parallel for private(i,j) collapse(2)
#endif
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            res[j+i*n] = v1[i]*v2[j];
        }
    }

    return;
}
void zero_array(double *v, int n) {
    int i;
    for(i=0;i<n;i++) {
        v[i] = 0;
    }
    return;
}


void cmatmat(double complex  *A, double complex *B, double complex *C,
					double complex alpha, double complex beta, int nA)
{
/* Performs \alpha * A.B + \beta * C and stores the output in C.
	A,B, and C are all matrices.
	This is essenitally a wrapper for the ZGEMM BLAS routine
*/
	int i,j;
	char TRANSA = 't';
	char TRANSB = 't';
	int m = nA;
	int n = nA;
	int k = nA;
	int LDA = nA;
	int LDB = nA;
	int LDC = nA;

    double complex *work = (double complex *)malloc(sizeof(double complex)*nA*nA);
	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++) work[i+nA*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = work[i+nA*j];
	}

	zgemm_(&TRANSA, &TRANSB, &m,&n,&k,&alpha,A,&LDA,B,&LDB,&beta,C,&LDC);


	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)  work[i+nA*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = work[i+nA*j];
	}
    SAFE_FREE(work);
	return;

}

void cmatvec(double complex  *A, double complex*B, double complex *C,
					double complex alpha, double complex beta, int nB)
{
/* Performs \alpha * A.B + \beta * C and stores the output in C.
	A is a matrix, B and C are vectors.
	This is essenitally a wrapper for the ZGEMV BLAS routine
*/

	char TRANS = 't';
	int m = nB;
	int n = nB;
	int LDA = nB;
	int INCX = 1;
	int INCY = 1;



	zgemv_(&TRANS, &m,&n,&alpha,A,&LDA,B,&INCX,&beta,C,&INCY);

	return;

}

void cmatmat3(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta) {

    C[0] = beta*C[0] + alpha*(A[0]*B[0] + A[1]*B[3] + A[2]*B[6]);
    C[1] = beta*C[1] + alpha*(A[0]*B[1] + A[1]*B[4] + A[2]*B[7]);
    C[2] = beta*C[2] + alpha*(A[0]*B[2] + A[1]*B[5] + A[2]*B[8]);

    C[3] = beta*C[3] + alpha*(A[3]*B[0] + A[4]*B[3] + A[5]*B[6]);
    C[4] = beta*C[4] + alpha*(A[3]*B[1] + A[4]*B[4] + A[5]*B[7]);
    C[5] = beta*C[5] + alpha*(A[3]*B[2] + A[4]*B[5] + A[5]*B[8]);

    C[6] = beta*C[6] + alpha*(A[6]*B[0] + A[7]*B[3] + A[8]*B[6]);
    C[7] = beta*C[7] + alpha*(A[6]*B[1] + A[7]*B[4] + A[8]*B[7]);
    C[8] = beta*C[8] + alpha*(A[6]*B[2] + A[7]*B[5] + A[8]*B[8]);

    return;
}

void cmatvec3(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta) {

    C[0] = beta*C[0] + alpha*(A[0]*B[0] + A[1]*B[1] + A[2]*B[2]);
    C[1] = beta*C[1] + alpha*(A[3]*B[0] + A[4]*B[1] + A[5]*B[2]);
    C[2] = beta*C[2] + alpha*(A[6]*B[0] + A[7]*B[1] + A[8]*B[2]);
    
    return;
}


void csolve(double complex *A, double complex *B,int nA, int nRHS) {

    int N = nA;
    int NRHS = nRHS;
    int LDA = nA;
    int *IPIV = (int *)malloc(sizeof(int)*nA);
    int LDB = nA;
    int INFO;

    double complex *AT = (double complex *)malloc(sizeof(double complex)*N*N);
    double complex *BT = (double complex *)malloc(sizeof(double complex)*N*nRHS);
    int i,j;
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            AT[i + N*j] = A[j + N*i];
        }
    }
    if (nRHS > 1) {
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++) {
                BT[i + N*j] = B[j + N*i];
            }
        }
    }
    else {
        for(i=0;i<N;i++) BT[i] = B[i];
    }

    zgesv_(&N, &NRHS, AT, &LDA, IPIV, BT, &LDB, &INFO);

    if (nRHS > 1) {
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++) {
                B[i + N*j] = BT[j + N*i];
            }
        }
    }
    else {
        for(i=0;i<N;i++) B[i] = BT[i];
    }
    SAFE_FREE(AT);
    SAFE_FREE(BT);
    return;
}

void cthomas_alg(double complex *a, double complex *b, double complex *c, double complex *d, int n) {
    /* 
     * Solve the tridiagonal system a_i*y_{i-1} + b_i*y_i + c_i*y_{i+1} = d_i
     * using Thomas's algorithm
     * Solution is stored in d array
     */
    int i; 

    c[0] /= b[0];
    d[0] /= b[0];

    for(i=1;i<n-1;i++) {
        
        b[i] -= a[i-1] * c[i-1];

        d[i] -= a[i-1]*d[i-1];


        c[i] = c[i] / b[i];
        d[i] = d[i] / b[i];
        

    }
    i = n-1;
    b[i] -= a[i-1] * c[i-1];
    d[i] -= a[i-1]*d[i-1];
    d[i] = d[i] / b[i];


    for(i=n-2;i>=0;i--) {
        d[i] -= c[i]*d[i+1];
    }

    return;
}
void cthomas_alg_block(double complex *a, double complex *b, double complex *c, double complex *d, int n, int m) {
    /* 
     * Solve the block tridiagonal system a_i*y_{i-1} + b_i*y_i + c_i*y_{i+1} = d_i
     * using Thomas's algorithm
     * Solution is stored in d array
     */
    int i; 
    int size = m*m;

    csolve(&b[0],&c[0],m,m);
    csolve(&b[0],&d[0],m,1);

    for(i=1;i<n-1;i++) {
        
        cmatmat3(&a[(i-1)*size],&c[(i-1)*size],&b[i*size],-1.0,1.0);
        cmatvec3(&a[(i-1)*size],&d[(i-1)*m],&d[i*m],-1.0,1.0);

        csolve(&b[i*size],&c[i*size],m,m);
        csolve(&b[i*size],&d[i*m],m,1);
        

    }
    i = n-1;
        
    cmatmat3(&a[(i-1)*size],&c[(i-1)*size],&b[i*size],-1.0,1.0);
    cmatvec3(&a[(i-1)*size],&d[(i-1)*m],&d[i*m],-1.0,1.0);
    csolve(&b[i*size],&d[i*m],m,1);


    for(i=n-2;i>=0;i--) {
        cmatvec3(&c[i*size],&d[(i+1)*m],&d[i*m],-1,1);
    }

    return;
}

void cconstruct_total_matrix(double complex *ld, double complex *md, double complex *ud, double complex *mat, int n, int m) {

    int i,j,k;

    int size = m*m;

    i=0;
    for(j=i*m;j<(i+1)*m;j++) {
        for(k=i*m; k < (i+1)*m; k++) {
            mat[k + j*n] = md[i*size + k + j*m];
        }
        for(k=(i+1)*m; k<(i+2)*m;k++) {
            mat[k + j*n] = ud[i*size + k + j*m];
        }
    }


    for(i=1;i<n-1;i++) {

        for(j=i*m;j<(i+1)*m;j++) {
            for(k=i*m; k < (i+1)*m; k++) {
                mat[k + j*n] = md[i*size + k + j*m];
            }
            for(k=(i+1)*m; k<(i+2)*m;k++) {
                mat[k + j*n] = ud[i*size + k + j*m];
            }
            for(k=(i-1)*m;k<i*m;k++) {
                mat[k + j*n] = ld[i*size + k+j*m];
            }
        }
    }
    i = n-1;
    for(j=i*m;j<(i+1)*m;j++) {
        for(k=i*m; k < (i+1)*m; k++) {
            mat[k + j*n] = md[i*size + k + j*m];
        }
        for(k=(i-1)*m;k<i*m;k++) {
            mat[k + j*n] = ld[i*size + k+j*m];
        }
    }


}


