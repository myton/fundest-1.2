/////////////////////////////////////////////////////////////////////////////////
// 
//  Fundamental matrix estimation from 2D point matches
//  Copyright (C) 2002 - 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "compiler.h"

#include "util.h"

/* lapack functions prototypes */

/* eigenvalues & eigenvectors (unpacked storage) */
extern int F77_FUNC(dsyev)(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);

/* SVD */
extern int F77_FUNC(dgesvd)(char *jobu, char *jobvt, int *m, int *n,
                            double *a, int *lda, double *s, double *u, int *ldu,
                            double *vt, int *ldvt, double *work, int *lwork, int *info);

/* lapack 3.0 routine, faster than dgesvd() */
extern int F77_FUNC(dgesdd)(char *jobz, int *m, int *n, double *a, int *lda,
                            double *s, double *u, int *ldu, double *vt, int *ldvt,
                            double *work, int *lwork, int *iwork, int *info);

/* matrix multiplication */
extern int F77_FUNC(dgemm)(char *transa, char *transb, int *m, int *n, int *k,
                          double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

/* LU decomposition, linear system solution and matrix inversion */
extern int F77_FUNC(dgetrf)(int *m, int *n, double *a, int *lda, int *ipiv, int *info); /* blocked LU */
extern int F77_FUNC(dgetf2)(int *m, int *n, double *a, int *lda, int *ipiv, int *info); /* unblocked LU */
extern int F77_FUNC(dgetri)(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

/* Solve min |Ax| subject to |x|=1
 * The solution is the eigenvector of A^T A corresponding
 * to the smallest eigenvalue.
 *
 * A is mxn, x is nx1
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int fundest_min_Ax_normEIG(double *A, int m, int n, double *x)
{
static double *buf=NULL;
static int buf_sz=0;

static int nb=0;

register int i, j, k;
double *eigvals, *eigvecs, *work, thresh;
int info, tot_sz, eigvals_sz, eigvecs_sz, work_sz;

  if(!A){
    if(buf) free(buf);
    buf=NULL;
    buf_sz=0;
    return 1;
  }

  /* calculate required memory size */
  if(!nb){
    /* determine optimal work size for dsyev */
    work_sz=-1; // workspace query; optimal size is returned in thresh
    F77_FUNC(dsyev)("V", "L", &n, NULL, &n, NULL, &thresh, &work_sz, &info);
    nb=((int)thresh)/n - 2; // optimal worksize is n*(nb+2);
  }
  eigvecs_sz=n*n;
  eigvals_sz=n;
  work_sz=(nb+2)*n;
  tot_sz=eigvecs_sz+eigvals_sz+work_sz;

  if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
    if(buf) free(buf); /* free previously allocated memory */

    buf_sz=tot_sz;
    buf=(double *)malloc(buf_sz*sizeof(double));
    if(!buf){
      fprintf(stderr, "Memory allocation request in fundest_min_Ax_normEIG() failed!\n");
      exit(1);
    }
  }

  eigvecs=buf;
  eigvals=eigvecs+eigvecs_sz;
  work=eigvals+eigvals_sz;

  /* initialize lower triangle of eigvecs to A^t * A */
#if 0
  for(j=0; j<n; j++){
    register double sum;

    for(i=j; i<n; i++){
      for(k=0, sum=0.0; k<m; k++)
        sum+=A[k*n+i]*A[k*n+j];
      eigvecs[j*n+i]=sum;
    }
  }
#else
  /* more cache-friendly scheme: accesses are always to consecutive addresses */
  memset(eigvecs, 0, n*n*sizeof(double));

  for(k=m; k-->0;  ){
    register double *Akn, alpha;

    Akn=A+k*n; // &A[k*n]
    for(j=n; j-->0;  ){
      alpha=Akn[j]; // A[k*n+j];
      for(i=j; i<n; i++)
        eigvecs[j*n+i]+=Akn[i]*alpha; // A[k*n+i]*alpha;
    }
  }
#endif

  F77_FUNC(dsyev)("V", "L", &n, (double *)eigvecs, &n, eigvals, work, &work_sz, &info);

  if(info<0){
    fprintf(stderr, "LAPACK error: illegal value for argument %d of dsyev in fundest_min_Ax_normEIG()\n", -info);
    exit(1);
  }
  else if(info>0){
    fprintf(stderr, "LAPACK error: dsyev failed to converge in fundest_min_Ax_normEIG();\n%d %s", info,
        "off-diagonal elements of an intermediate tridiagonal form did not converge to zero\n");

    return 0;
  }

  for(i=n-1, j=0, thresh=eigvals[n-1]*DBL_EPSILON; i>0; --i, ++j)
    if(eigvals[i]<=thresh) break; /* remaining eigenvalues are smaller than this */


  if(j!=n-1){ /* matrix is not of rank n-1! */
    //fprintf(stderr, "Unacceptable rank %d in fundest_min_Ax_normEIG\n", rank);

    return 0;
  }

  /* min eigenvalue is the first (they are returned in ascending order), return the first eigenvector */
  for(i=0; i<n; i++)
    x[i]=eigvecs[i];

  return 1;
}

/* Solve min |Ax| subject to |x|=1
 * The solution is the right singular vector (Vt) corresponding to A's
 * minimum singular value: A=U*D*V^T ==> A^T*A=V*D^2*V^T
 *
 * A is mxn, x is nx1
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int fundest_min_Ax_normSVD(double *A, int m, int n, double *x)
{
static double *buf=NULL;
static int buf_sz=0;

register int i, j;
double *a, *u, *s, *vt, *work, thresh;
int info, worksz, *iwork, iworksz;
int a_sz, u_sz, s_sz, vt_sz, tot_sz;
//static double eps=-1.0;

  if(!A){
    if(buf) free(buf);
    buf=NULL;
    buf_sz=0;
    return 1;
  }

  /* calculate required memory size. Note that if dgesvd is used, the memory for u is not actually needed... */
  worksz=-1; // workspace query. Keep in mind that dgesdd requires more memory than dgesvd
  /* note that optimal work size is returned in thresh */
  //F77_FUNC(dgesdd)("A", (int *)&m, (int *)&n, NULL, (int *)&m, NULL, NULL, (int *)&m, NULL, (int *)&n, (double *)&thresh, (int *)&worksz, NULL, &info);
  F77_FUNC(dgesvd)("N", "A", (int *)&m, (int *)&n, NULL, (int *)&m, NULL, NULL, (int *)&m, NULL, (int *)&n, (double *)&thresh, (int *)&worksz, &info);
  worksz=(int)thresh;
  iworksz=8*n;
  a_sz=m*n;
  u_sz=m*m; s_sz=n; vt_sz=n*n;

  tot_sz=(a_sz + u_sz + s_sz + vt_sz + worksz)*sizeof(double) + iworksz*sizeof(int); /* should be arranged in that order for proper doubles alignment */

  if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
    if(buf) free(buf); /* free previously allocated memory */

    buf_sz=tot_sz;
    buf=(double *)malloc(buf_sz);
    if(!buf){
      fprintf(stderr, "Memory allocation request for %d bytes in fundest_min_Ax_normSVD() failed!\n", buf_sz);
      exit(1);
    }
  }

  a=buf;
  u=a+a_sz;
  s=u+u_sz;
  vt=s+s_sz;
  work=vt+vt_sz;
  iwork=(int *)(work+worksz);

  /* store A (column major!) into a */
  for(i=0; i<m; i++)
    for(j=0; j<n; j++)
      a[i+j*m]=A[i*n+j];

  /* SVD decomposition of A */
  //F77_FUNC(dgesdd)("A", (int *)&m, (int *)&n, a, (int *)&m, s, u, (int *)&m, vt, (int *)&n, work, (int *)&worksz, iwork, &info);
  F77_FUNC(dgesvd)("N", "A", (int *)&m, (int *)&n, a, (int *)&m, s, NULL, (int *)&m, vt, (int *)&n, work, (int *)&worksz, &info);

  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dgesvd/dgesdd in fundest_min_Ax_normSVD()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "LAPACK error: dgesdd (dbdsdc)/dgesvd (dbdsqr) failed to converge in fundest_min_Ax_normSVD() [info=%d]\n", info);
      return 0;
    }
  }

#if 0
  /* compute machine epsilon */
  if(eps<0.0){
    double aux;

    for(eps=1.0; aux=eps+1.0, aux-1.0>0.0; eps*=0.5)
                              ;
    eps*=2.0;
  }
#endif

  /* determine A's rank */
  j=(m<=n)? m : n; // min(m, n);
  for(i=0, thresh=DBL_EPSILON*s[0]; i<j; ++i)
    if(s[i]<=thresh) break; /* remaining singular values are smaller than this */

  if(i<n-1){
    //fprintf(stderr, "Unacceptable rank %d in fundest_min_Ax_normSVD()\n", i);
    memset(x, 0, n*sizeof(double));
    return 0; /* A should have rank n-1 */
  }

  /* s[n-1] is the smallest singular value */
  vt+=n-1;
  for(j=0; j<n; ++j)
    x[j]=vt[j*n]; //vt[n-1+j*n];

  return 1;
}

/* Solve min |Ax| subject to |x|=1 and x=G*x' for a given rank r matrix G and an unknown vector x',
 * i.e. min |AGx'| subject to |Gx'|=1
 *
 * The solution is as described in A5.4.1 of HZ2 (p.595)
 *
 * A is mxn, x is nx1 and G nxn
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int fundest_min_AGx_norm(double *A, double *G, int m, int n, int r, double *x)
{
static double *buf=NULL;
static int buf_sz=0;

register int i, j;
double *g, *u, *s, *vt, *work, thresh;
double *Aup, *y, *up;
int info, worksz, *iwork, iworksz;
int g_sz, u_sz, s_sz, vt_sz, Aup_sz, y_sz, tot_sz;
//static double eps=-1.0;
double one=1.0, zero=0.0;

  if(!A){
    if(buf) free(buf);
    buf=NULL;
    buf_sz=0;
    return 1;
  }
  r=(r<=n)? r : n; // min(r, n);

  /* calculate required memory size */
  worksz=-1; // workspace query
  /* note that optimal work size is returned in thresh */
  //F77_FUNC(dgesdd)("A", (int *)&n, (int *)&n, NULL, (int *)&n, NULL, NULL, (int *)&n, NULL, (int *)&n, (double *)&thresh, (int *)&worksz, NULL, &info);
  F77_FUNC(dgesvd)("A", "N", (int *)&n, (int *)&n, NULL, (int *)&n, NULL, NULL, (int *)&n, NULL, (int *)&n, (double *)&thresh, (int *)&worksz, &info);
  worksz=(int)thresh;
  iworksz=8*n;
  g_sz=n*n;
  u_sz=n*n; s_sz=n; vt_sz=n*n;
  if(r>0){
    Aup_sz=m*r;
    y_sz=r;
  }
  else{ // r not yet known but n>=r
    Aup_sz=m*n;
    y_sz=n;
  }

  tot_sz=(g_sz + u_sz + s_sz + vt_sz + Aup_sz + y_sz + worksz)*sizeof(double) + iworksz*sizeof(int); /* should be arranged in that order for proper doubles alignment */

  if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
    if(buf) free(buf); /* free previously allocated memory */

    buf_sz=tot_sz;
    buf=(double *)malloc(buf_sz);
    if(!buf){
      fprintf(stderr, "Memory allocation request for %d bytes in fundest_min_AGx_norm() failed!\n", buf_sz);
      exit(1);
    }
  }

  g=buf;
  u=g+g_sz;
  s=u+u_sz;
  vt=s+s_sz;
  Aup=vt+vt_sz;
  y=Aup+Aup_sz;
  work=y+y_sz;
  iwork=(int *)(work+worksz);

  /* store G (column major!) into g */
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      g[i+j*n]=G[i*n+j];

  /* SVD decomposition of G */
  //F77_FUNC(dgesdd)("A", (int *)&n, (int *)&n, g, (int *)&n, s, u, (int *)&n, vt, (int *)&n, work, (int *)&worksz, iwork, &info);
  F77_FUNC(dgesvd)("A", "N", (int *)&n, (int *)&n, g, (int *)&n, s, u, (int *)&n, vt, (int *)&n, work, (int *)&worksz, &info);

  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dgesvd/dgesdd in fundest_min_AGx_norm()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "LAPACK error: dgesdd (dbdsdc)/dgesvd (dbdsqr) failed to converge in fundest_min_AGx_norm() [info=%d]\n", info);
      return 0;
    }
  }

#if 0
  /* compute machine epsilon */
  if(eps<0.0){
    double aux;

    for(eps=1.0; aux=eps+1.0, aux-1.0>0.0; eps*=0.5)
                              ;
    eps*=2.0;
  }
#endif

  if(r<=0){ /* determine G's rank */
    for(i=0, thresh=DBL_EPSILON*s[0]; i<n; ++i)
      if(s[i]<=thresh) break; /* remaining singular values are smaller than this */
    r=i;
  }

  /* keep the r first columns of u in up and clear the rest */
  //for(i=r*n; i<n*n; ++i) u[i]=0.0;
  memset(u+r*n, 0, n*(n-r)*sizeof(double));
  up=u;

  /* compute the mxr product A*up in *row major*. Note that A is in row major
   * whereas up in column major. dgemm() is asked to compute up'*A' (i.e the product in row maj.),
   * so it has to transpose up and use A as is (i.e. transpose in col. maj.)
   */
  F77_FUNC(dgemm)("T", "N", (int *)&r, (int *)&m, (int *)&n, (double *)&one, up, (int *)&n, A, (int *)&n, (double *)&zero, Aup, (int *)&r);

  /* solve min |A*up| with |y|=1 */
  fundest_min_Ax_normEIG(Aup, m, r, y);

  /* x=up*y */
  for(i=0; i<n; ++i){
    register double sum;

    for(j=0, sum=0.0; j<r; ++j)
      sum+=up[j*n+i]*y[j];
    x[i]=sum;
  }

  /* x' is v'*d'^-1*y; note that to compute V, the 2nd argument of dgesvd should be changed to "A"! */

  return 1;
}

/*
 * This function computes the inverse of a square matrix A into B
 * using LU decomposition
 *
 * The function returns 0 in case of error (e.g. A is singular),
 * 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int fundest_matinvLU(double *A, double *B, int m)
{
static double *buf=NULL;
static int buf_sz=0, nb=0;

int a_sz, ipiv_sz, work_sz, tot_sz;
register int i, j;
int info, *ipiv;
double *a, *work;
   
  if(A==NULL){
    if(buf) free(buf);
    buf=NULL;
    buf_sz=0;

    return 1;
  }

  /* calculate required memory size */
  ipiv_sz=m;
  a_sz=m*m;
  if(!nb){
    double tmp;

    work_sz=-1; // workspace query; optimal size is returned in tmp
    F77_FUNC(dgetri)((int *)&m, NULL, (int *)&m, NULL, (double *)&tmp, (int *)&work_sz, (int *)&info);
    nb=((int)tmp)/m; // optimal worksize is m*nb
  }
  work_sz=nb*m;
  tot_sz=(a_sz + work_sz)*sizeof(double) + ipiv_sz*sizeof(int); /* should be arranged in that order for proper doubles alignment */

  if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
    if(buf) free(buf); /* free previously allocated memory */

    buf_sz=tot_sz;
    buf=(double *)malloc(buf_sz);
    if(!buf){
      fprintf(stderr, "memory allocation in fundest_matinvLU() failed!\n");
      exit(1);
    }
  }

  a=buf;
  work=a+a_sz;
  ipiv=(int *)(work+work_sz);

  /* store A (column major!) into a */
	for(i=0; i<m; ++i)
		for(j=0; j<m; ++j)
			a[i+j*m]=A[i*m+j];

  /* LU decomposition for A */
	//F77_FUNC(dgetrf)((int *)&m, (int *)&m, a, (int *)&m, ipiv, (int *)&info);  
	F77_FUNC(dgetf2)((int *)&m, (int *)&m, a, (int *)&m, ipiv, (int *)&info);  
	if(info!=0){
		if(info<0){
			fprintf(stderr, "argument %d of dgetf2/dgetrf illegal in fundest_matinvLU()\n", -info);
			exit(1);
		}
		else{
			fprintf(stderr, "singular matrix A for dgetf2/dgetrf in fundest_matinvLU()\n");
			return 0;
		}
	}

  /* (A)^{-1} from LU */
	F77_FUNC(dgetri)((int *)&m, a, (int *)&m, ipiv, work, (int *)&work_sz, (int *)&info);
	if(info!=0){
		if(info<0){
			fprintf(stderr, "argument %d of dgetri illegal in fundest_matinvLU()\n", -info);
			exit(1);
		}
		else{
			fprintf(stderr, "singular matrix A for dgetri in fundest_matinvLU()\n");
			return 0;
		}
	}

	/* store (A)^{-1} in B */
	for(i=0; i<m; ++i)
		for(j=0; j<m; ++j)
      B[i*m+j]=a[i+j*m];

	return 1;
}
