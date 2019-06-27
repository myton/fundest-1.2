/////////////////////////////////////////////////////////////////////////////////
//
//  Fundamental matrix estimation from 2D point matches
//  Copyright (C) 2002 - 2015  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  An initial estimate is computed with the normalized 8 point algorithm and
//  then refined nonlinearly by minimizing the symmetric distance to epipolar
//  lines using a parameterization which enforces the rank 2 constraint
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
#include "fundest.h"
#include "calc_fund_coeffs.h"

#include "lqs.h"
//#include "ransac.h"

#include <levmar.h>

//#define USE_BUCKETS // CHECKME


/* minimal number of matches for the linear algorithm */
#define NUM_LINFMATCHES    8 // 8-point algorithm

/* number of parameters for the nonlinear refinement */
#define NUM_FNLPARAMS 8

#define SQR(x)      ((x)*(x))
#define FABS(x)     (((x)>=0)? (x) : -(x))

/* parameterizations for F: rows _xy make up the remaining row */
#define _FPARAM_12    12
#define _FPARAM_13    13
#define _FPARAM_23    23

/* normalization */
#define _DONORM   1

/* SVD */
extern int F77_FUNC(dgesvd)(
        char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s,
        double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);



/***** robust estimation of epipolar geometry */

/* variables used by various estimation routines */
struct Fdata{
  double **cfs;
  double (*pts0)[2], (*pts1)[2];
  int *inliersidx, numInliers;

  /* following elements are used for avoiding multiple mallocs in estLinEpipGeom() */
  double *M;
  int Mrows;
};

/* unrolled memory copy for 9-sized arrays */
inline static void dmemcpy9(double *dest, double *src)
{
  *dest++=*src++;
  *dest++=*src++;
  *dest++=*src++;
  *dest++=*src++;
  *dest++=*src++;
  *dest++=*src++;
  *dest++=*src++;
  *dest++=*src++;
  *dest=*src;
}

/* linearly estimate F from point matches s.t. m_1^t * F * m_2=0, with m_0, m_1 specified by ptsidx.
 * if ptsidx==NULL, all points are employed
 */
static int estLinEpipGeom(double *p, int npts, int *ptsidx, void *adata)
{
register int i;
double *mat;
register double *matrow;
static int (*const solver)(double *, int, int, double *)=
                    fundest_min_Ax_normSVD; /* using SVD is more stable but slower & needs more memory */
                    //fundest_min_Ax_normEIG;
struct Fdata *dat=(struct Fdata *)adata;

  if(p==NULL){ /* cleanup */
    (*solver)(NULL, 0, 0, NULL);

    /* free working memory */
    free(dat->M);
    dat->M=NULL; dat->Mrows=0;
    return 0;
  }

  if(dat->Mrows!=npts){ /* free existing matrix and allocate a new one */
    if(dat->M) /* is there anything to free? */
      free(dat->M);

    if((mat=(double *)malloc(npts*NUM_FPARAMS*sizeof(double)))==NULL){
      fprintf(stderr, "Memory allocation request failed in estLinEpipGeom()\n");
      exit(1);
    }
    dat->M=mat;
    dat->Mrows=npts;
  }

  mat=dat->M;

  /* fill in constraints matrix row by row */
  if(ptsidx){
    for(i=0, matrow=mat; i<npts; ++i, matrow+=NUM_FPARAMS)
      dmemcpy9(matrow, dat->cfs[ptsidx[i]]);
  }
  else{
    for(i=0, matrow=mat; i<npts; ++i, matrow+=NUM_FPARAMS)
      dmemcpy9(matrow, dat->cfs[i]);
  }

  /* solve min |mat*p|, |p|=1 */
  if(!(*solver)(mat, npts, NUM_FPARAMS, p))
    return 0;

  return 1;
}

/* compute the (squared) algebraic residuals corresponding to F */
static void fmatrixLinResidualsAlg(double *F, int numres, void *adata, double *resid)
{
register int i, j;
double *cf, sum;
struct Fdata *dat=(struct Fdata *)adata;

  for(i=0; i<numres; ++i){
    cf=dat->cfs[i];
    sum=0.0;
    for(j=0; j<NUM_FPARAMS; ++j){
      sum+=cf[j]*F[j];
    }
    resid[i]=sum*sum;
  }
}

/* compute the geometric residuals corresponding to F as the (squared) symmetric
 * distance between a point and the epipolar line
 */
static void fmatrixLinResidualsGeom(double *F, int numres, void *adata, double *resid)
{
register int i;
double *pt0, *pt1, line[3], tmp;
struct Fdata *dat=(struct Fdata *)adata;

  for(i=0; i<numres; ++i){
    pt0=dat->pts0[i]; pt1=dat->pts1[i];

    /* img 1 */
    line[0]=F[0]*pt0[0]+ F[1]*pt0[1]+ F[2];
    line[1]=F[3]*pt0[0]+ F[4]*pt0[1]+ F[5];
    line[2]=F[6]*pt0[0]+ F[7]*pt0[1]+ F[8];

    tmp=(line[0]*pt1[0] + line[1]*pt1[1] + line[2]);
    resid[i]=SQR(tmp)/(SQR(line[0]) + SQR(line[1]));

    /* img 0 */
    line[0]=F[0]*pt1[0]+ F[3]*pt1[1]+ F[6];
    line[1]=F[1]*pt1[0]+ F[4]*pt1[1]+ F[7];
    line[2]=F[2]*pt1[0]+ F[5]*pt1[1]+ F[8];
    tmp=(line[0]*pt0[0] + line[1]*pt0[1] + line[2]);
    resid[i]+=SQR(tmp)/(SQR(line[0]) + SQR(line[1]));
  }
}


/** Distance from epipolar lines **/

/* epipolar line error and Jacobian (12) */
static void epipDist_12(double *p, double *x, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;

  for(i=0; i<ninl; ++i){
    calcEpipDist_12(pts0[inlidx[i]], pts1[inlidx[i]], p, x+i);
  }
}

static void epipDistJac_12(double *p, double *jac, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;
register double *jacrow;

  //memset(jac, 0, m*n*sizeof(double));

  for(i=0, jacrow=jac; i<ninl; ++i, jacrow+=m){
    calcEpipDistJac_12(pts0[inlidx[i]], pts1[inlidx[i]], p, jacrow);
  }
}

/* epipolar line error and Jacobian (13) */
static void epipDist_13(double *p, double *x, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;

  for(i=0; i<ninl; ++i){
    calcEpipDist_13(pts0[inlidx[i]], pts1[inlidx[i]], p, x+i);
  }
}

static void epipDistJac_13(double *p, double *jac, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;
register double *jacrow;

  //memset(jac, 0, m*n*sizeof(double));

  for(i=0, jacrow=jac; i<ninl; ++i, jacrow+=m){
    calcEpipDistJac_13(pts0[inlidx[i]], pts1[inlidx[i]], p, jacrow);
  }
}

/* epipolar line error and Jacobian (23) */
static void epipDist_23(double *p, double *x, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;

  for(i=0; i<ninl; ++i){
    calcEpipDist_23(pts0[inlidx[i]], pts1[inlidx[i]], p, x+i);
  }
}

static void epipDistJac_23(double *p, double *jac, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;
register double *jacrow;

  //memset(jac, 0, m*n*sizeof(double));

  for(i=0, jacrow=jac; i<ninl; ++i, jacrow+=m){
    calcEpipDistJac_23(pts0[inlidx[i]], pts1[inlidx[i]], p, jacrow);
  }
}

/** Sampson **/

/* Sampson error and Jacobian (12) */
static void epipConstrSampson_12(double *p, double *x, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;

  for(i=0; i<ninl; ++i){
    calcEpipConstrSampsonDist_12(pts0[inlidx[i]], pts1[inlidx[i]], p, x+i);
  }
}

static void epipConstrSampsonJac_12(double *p, double *jac, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;
register double *jacrow;

  //memset(jac, 0, m*n*sizeof(double));

  for(i=0, jacrow=jac; i<ninl; ++i, jacrow+=m){
    calcEpipConstrSampsonDistGrads_12(pts0[inlidx[i]], pts1[inlidx[i]], p, jacrow);
  }
}

/* Sampson error and Jacobian (13) */
static void epipConstrSampson_13(double *p, double *x, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;

  for(i=0; i<ninl; ++i){
    calcEpipConstrSampsonDist_13(pts0[inlidx[i]], pts1[inlidx[i]], p, x+i);
  }
}

static void epipConstrSampsonJac_13(double *p, double *jac, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;
register double *jacrow;

  //memset(jac, 0, m*n*sizeof(double));

  for(i=0, jacrow=jac; i<ninl; ++i, jacrow+=m){
    calcEpipConstrSampsonDistGrads_13(pts0[inlidx[i]], pts1[inlidx[i]], p, jacrow);
  }
}

/* Sampson error and Jacobian (23) */
static void epipConstrSampson_23(double *p, double *x, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;

  for(i=0; i<ninl; ++i){
    calcEpipConstrSampsonDist_23(pts0[inlidx[i]], pts1[inlidx[i]], p, x+i);
  }
}

static void epipConstrSampsonJac_23(double *p, double *jac, int m, int n, void *adata)
{
register int i;
struct Fdata *dat=(struct Fdata *)adata;
int ninl=dat->numInliers, *inlidx=dat->inliersidx;
double (*pts0)[2]=dat->pts0, (*pts1)[2]=dat->pts1;
register double *jacrow;

  //memset(jac, 0, m*n*sizeof(double));

  for(i=0, jacrow=jac; i<ninl; ++i, jacrow+=m){
    calcEpipConstrSampsonDistGrads_23(pts0[inlidx[i]], pts1[inlidx[i]], p, jacrow);
  }
}


/* determine r, s s.t. v2=r*v0 + s*v1.
 * The solution is computed in the least squares sense, using the
 * normal equations. Returns the absolute value of the determinant of A^T*A,
 * where A:=array(1..3, 1..2, [[v0[1], v1[1]], [v0[2], v1[2]], [v0[3], v1[3]]]);
 */
static double normeq3x3(double v0[3], double v1[3], double v2[3], double *r, double *s)
{
double a1, a2, b1, b2, c1, c2;
double det, fdet;

  /* Cramer's rule:
   * if  a1*r + b1*s=c1
   *     a2*r + b2*s=c2
   * then r=((c1*b2)-(c2*b1))/((a1*b2)-(a2*b1))
   *      s=((a1*c2)-(a2*c1))/((a1*b2)-(a2*b1))
   */
  /* A^T*A */
  a1=   v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2];
  b1=a2=v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
  b2=   v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];

  /* A^T*v2 */
  c1=v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2];
  c2=v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

#if 0
// use only two of the three available equations
  a1=v0[0]; b1=v1[0]; c1=v2[0]; // column 0
  a2=v0[1]; b2=v1[1]; c2=v2[1]; // column 1
  //a2=v0[2]; b2=v1[2]; c2=v2[2]; // column 2
#endif

  det=(a1*b2)-(a2*b1);
  fdet=FABS(det);

  if(fdet<1E-22){
    *r=*s=0.0;
    return 0.0;
  }

  det=1.0/det;
  *r=((c1*b2)-(c2*b1))*det;
  *s=((a1*c2)-(a2*c1))*det;

  return fdet;
}

/* unrolled memory copy functions for arrays sized 3 & 6 */
inline static void dmemcpy3(double *dest, double *src)
{
  *dest++=*src++;
  *dest++=*src++;
  *dest=*src;
}

inline static void dmemcpy6(double *dest, double *src)
{
  *dest++=*src++;
  *dest++=*src++;
  *dest++=*src++;
  *dest++=*src++;
  *dest++=*src++;
  *dest=*src;
}

/* nonlinear refinement of a fundamental matrix as specified by "howto":
 *  FUNDEST_EPIP_DIST: minimize the symmetric epipolar line distance
 *  FUNDEST_SAMPSON_ERROR: minimize the Sampson distance
 */
static void refineEpipGeom(double F[NUM_FPARAMS], int howto, struct Fdata *data, int verbose)
{
double opts[LM_OPTS_SZ], info[LM_INFO_SZ], f8[NUM_FNLPARAMS];
int m, n; // # unknowns & # constraints
int ninl=data->numInliers;
void (*err)(double *p, double *hx, int m, int n, void *adata)=NULL;
void (*jacerr)(double *p, double *j, int m, int n, void *adata)=NULL;
double fdet_12, fdet_13, fdet_23, r_12, s_12, r_13, s_13, r_23, s_23;
int ptyp, ret;

static char *howtoname[]={"FUNDEST_NO_EST", "FUNDEST_8PT", "FUNDEST_ALGMIN", "FUNDEST_EPIP_DIST", "FUNDEST_SAMPSON_ERROR"};

  opts[0]=LM_INIT_MU; opts[1]=1E-12; opts[2]=1E-12; opts[3]=1E-15;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  m=NUM_FNLPARAMS; n=ninl;

  /* determine F's parameterization so that one of its rows is a
   * linear combination of the other two:
   */
    /* r & s s.t. r*F(1, :) + s*F(2, :) = F(3, :) */
  fdet_12=normeq3x3(F,   F+3, F+6, &r_12, &s_12);
    /* r & s s.t. r*F(1, :) + s*F(3, :) = F(2, :) */
  fdet_13=normeq3x3(F,   F+6, F+3, &r_13, &s_13);
    /* r & s s.t. r*F(2, :) + s*F(3, :) = F(1, :) */
  fdet_23=normeq3x3(F+3, F+6, F,   &r_23, &s_23);

  if(fdet_12+fdet_13+fdet_23==0.0){
    fprintf(stderr, "*** cannot express any row of F as a linear combination of the other two in refineEpipGeom() -- zero matrix?\n");
    fprintf(stderr, "%.8lf %.8lf %.8lf\n", F[0], F[1], F[2]);
    fprintf(stderr, "%.8lf %.8lf %.8lf\n", F[3], F[4], F[5]);
    fprintf(stderr, "%.8lf %.8lf %.8lf\n", F[6], F[7], F[8]);

    return;
  }

  /* select the parameterization giving rise to the largest |determinant| */
  if(fdet_12>=fdet_13){
    if(fdet_12>=fdet_23){ /* 12 */
      dmemcpy6(f8, F); // copy the first two rows
      f8[6]=r_12; f8[7]=s_12; // add the two coefficients
      ptyp=_FPARAM_12;
    }else{ /* 23 */
      dmemcpy6(f8, F+3); // copy the last two rows
      f8[6]=r_23; f8[7]=s_23; // add the two coefficients
      ptyp=_FPARAM_23;
    }
  }else{
    if(fdet_13>=fdet_23){ /* 13 */
      dmemcpy3(f8  , F); // copy first row
      dmemcpy3(f8+3, F+6); // copy third row
      f8[6]=r_13; f8[7]=s_13; // add the two coefficients
      ptyp=_FPARAM_13;
    }else{ /* 23 */
      dmemcpy6(f8, F+3); // copy the last two rows
      f8[6]=r_23; f8[7]=s_23; // add the two coefficients
      ptyp=_FPARAM_23;
    }
  }

  //ptyp=_FPARAM_12; // enforce a certain parameterization, e.g. 12

  if(verbose){
    fprintf(stdout, "Parameterization chosen for refinement: %d (dets: %g %g %g)\n", ptyp, fdet_12, fdet_13, fdet_23);
  }

  if(howto==FUNDEST_SAMPSON_ERROR){
    switch(ptyp){
      case _FPARAM_12:
        err=epipConstrSampson_12;
        jacerr=epipConstrSampsonJac_12;
      break;

      case _FPARAM_13:
        err=epipConstrSampson_13;
        jacerr=epipConstrSampsonJac_13;
      break;

      case _FPARAM_23:
        err=epipConstrSampson_23;
        jacerr=epipConstrSampsonJac_23;
      break;

      default: /* should not reach this point */
        fprintf(stderr, "Intermal error #2 in refineEpipGeom()!\n");
        exit(1);
    }
  }
  else{ // ptyp==FUNDEST_EPIP_DIST
    switch(ptyp){
      case _FPARAM_12:
        err=epipDist_12;
        jacerr=epipDistJac_12;
      break;

      case _FPARAM_13:
        err=epipDist_13;
        jacerr=epipDistJac_13;
      break;

      case _FPARAM_23:
        err=epipDist_23;
        jacerr=epipDistJac_23;
      break;

      default: /* should not reach this point */
        fprintf(stderr, "Intermal error #1 in refineEpipGeom()!\n");
        exit(1);
    }
  }

  /* note that measurement vector x is NULL (i.e. 0) below */
#if 1
  ret=dlevmar_der(err, jacerr, f8, NULL, m, n, 1000, opts, info, NULL, NULL, (void *)data); // with analytic Jacobian
  //ret=dlevmar_dif(err, f8, NULL, m, n, 1000, opts, info, NULL, NULL, (void *)data); // no Jacobian
#else
  { double rp[2];

  rp[0]=LM_HUBER; rp[1]=0.5; // Huber
  //rp[0]=LM_TUKEY; rp[1]=0.005; // Tukey
  ret=dlevmar_rob_der(err, jacerr, f8, NULL, m, n, rp, 1000, opts, info, NULL, NULL, (void *)data); // with analytic Jacobian
  }
#endif

  switch(ptyp){
    double r, s;

  case _FPARAM_12:
    dmemcpy6(F, f8); // copy back the first two rows
    /* set F(3, :) = r*F(1, :) + s*F(2, :) */
    r=f8[6]; s=f8[7];
    F[6]=r*F[0] + s*F[3];
    F[7]=r*F[1] + s*F[4];
    F[8]=r*F[2] + s*F[5];
    break;

  case _FPARAM_13:
    dmemcpy3(F,   f8); // copy back row 1
    dmemcpy3(F+6, f8+3); // copy back row 3
    /* set F(2, :) = r*F(1, :) + s*F(3, :) */
    r=f8[6]; s=f8[7];
    F[3]=r*F[0] + s*F[6];
    F[4]=r*F[1] + s*F[7];
    F[5]=r*F[2] + s*F[8];
    break;

  case _FPARAM_23:
    dmemcpy6(F+3, f8); // copy back the last two rows
    /* set F(1, :) = r*F(2, :) + s*F(3, :) */
    r=f8[6]; s=f8[7];
    F[0]=r*F[3] + s*F[6];
    F[1]=r*F[4] + s*F[7];
    F[2]=r*F[5] + s*F[8];
    break;

  default: /* should not reach this point */
    fprintf(stderr, "Intermal error #3 in refineEpipGeom()!\n");
    exit(1);
  }

  if(verbose){
    fprintf(stdout, "\nRefinement method %s, %d measurements, %d variables\n", howtoname[howto], n, m);
    fprintf(stdout, "LM returned %d in %g iter, reason %g, error %g [initial %g], %d/%d func/fjac evals\n",
                    ret, info[5], info[6], info[1]/ninl, info[0]/ninl, (int)info[7], (int)info[8]);
#if 0
    fprintf(stdout, "\nSolution: ");
    for(i=0; i<m; ++i)
      fprintf(stdout, "%.7g ", f8[i]);
    fprintf(stdout, "\n");
    fprintf(stdout, "Parameterization: %d\n", ptyp);
#endif
  }
}

/* force F to have rank 2 */
static void mk2Rank3x3(double F[9])
{
register int i, j;
double svals[3], Umat[9], Vtmat[9], sum;

  /* transpose F in place to make it column major */
  sum=F[1]; F[1]=F[3]; F[3]=sum;
  sum=F[5]; F[5]=F[7]; F[7]=sum;
  sum=F[6]; F[6]=F[2]; F[2]=sum;

#if 1 // generic SVD
  {
  const int three=3;
  double work[32];
  const int worksz=32; /* more than needed */
  int info;

  F77_FUNC(dgesvd)("A", "A", (int*)&three, (int*)&three, F, (int*)&three, svals, Umat, (int*)&three,
                   Vtmat, (int*)&three, work, (int*)&worksz, &info);
  if(info<0){
    fprintf(stderr, "LAPACK error: illegal value for argument %d of dgesvd in mk2Rank3x3()\n", -info);
    exit(1);
  }
  else if(info>0){
    fprintf(stderr, "LAPACK error: dgesvd (dbdsqr) failed to converge in mk2Rank3x3();\n%d %s", info,
                    "superdiagonals of an intermediate bidiagonal form B did not converge to zero\n");
    exit(1);
  }
  }
#else // custom 3x3 SVD
  svd3(Umat, svals, Vtmat, F);
  /* svd3 actually returns V; for compatibility with LAPACK which returns V^T,
   * the computed Vt is transposed in place to yield the true Vt
   */
  sum=Vtmat[1]; Vtmat[1]=Vtmat[3]; Vtmat[3]=sum;
  sum=Vtmat[5]; Vtmat[5]=Vtmat[7]; Vtmat[7]=sum;
  sum=Vtmat[6]; Vtmat[6]=Vtmat[2]; Vtmat[2]=sum;
#endif

  svals[2]=0.0; /* singular values have been sorted in descending order */
  /* compute the rank-2 matrix */
  for(i=0; i<3; i++)
    for(j=0; j<3; j++){
      sum =Umat[i+0*3]*svals[0]*Vtmat[j*3+0];
      sum+=Umat[i+1*3]*svals[1]*Vtmat[j*3+1];
      sum+=Umat[i+2*3]*svals[2]*Vtmat[j*3+2];
      F[i*3+j]=sum;
    }
}

/* normalize 3x3 matrix M so that its infinity norm is 1 */
static void inftynorm3x3(double M[9])
{
double sum, max;

  max=FABS(M[0])+FABS(M[1])+FABS(M[2]);

  sum=FABS(M[3])+FABS(M[4])+FABS(M[5]);
  if(sum>max) max=sum;

  sum=FABS(M[6])+FABS(M[7])+FABS(M[8]);
  if(sum>max) max=sum;

  if(max<1E-12){
    fprintf(stderr, "Zero matrix in inftynorm3x3()\n");
    return;
  }

  max=1.0/max;

  M[0]*=max; M[1]*=max; M[2]*=max;
  M[3]*=max; M[4]*=max; M[5]*=max;
  M[6]*=max; M[7]*=max; M[8]*=max;
}

/*********************************************** minimization of algebraic error ***********************************************/
struct Fdatalg{
  double prevepip[3]; // epipole from previous iteration
  double fm[NUM_FPARAMS]; // storage for the estimated F
  double **cfs;
};

/* compute the algebraic error, i.e. steps (ii) & (iii) of Alg. 11.2 */
static void algerror(double *e, double *x, int m, int n, void *adata)
{
register int i;
struct Fdatalg *datalg=(struct Fdatalg *)adata;
register double **A=datalg->cfs, *f=datalg->fm;
double *ep=datalg->prevepip, E[9*9];
int sg;

  sg=(e[0]*ep[0]+e[1]*ep[1]+e[2]*ep[2]<0.0)? -1 : 1; // sign check

  /* setup E=diag([e]x [e]x [e]x) */
  memset(E, 0, 9*9*sizeof(double));

  E[1]=-e[2]; E[2]=e[1];
  E[9]=e[2]; E[11]=-e[0];
  E[18]=-e[1]; E[19]=e[0];

  E[1+30]=-e[2]; E[2+30]=e[1];
  E[9+30]=e[2]; E[11+30]=-e[0];
  E[18+30]=-e[1]; E[19+30]=e[0];

  E[1+60]=-e[2]; E[2+60]=e[1];
  E[9+60]=e[2]; E[11+60]=-e[0];
  E[18+60]=-e[1]; E[19+60]=e[0];

  fundest_min_AGx_norm(A[0], E, n, NUM_FPARAMS, 6, f);

  // flip sign
  if(sg<0){
    f[0]=-f[0]; f[1]=-f[1]; f[2]=-f[2]; f[3]=-f[3];
    f[4]=-f[4]; f[5]=-f[5]; f[6]=-f[6]; f[7]=-f[7]; f[8]=-f[8];
  }

  /* x=A*f */
  for(i=0; i<n; ++i){
    double *a=A[i];

    x[i]=(a[0]*f[0] + a[1]*f[1]) + (a[2]*f[2] + a[3]*f[3]) +
         (a[4]*f[4] + a[5]*f[5]) + (a[6]*f[6] + a[7]*f[7]) + a[8]*f[8];
  }
 
  // prepare for next iteration
  ep[0]=e[0]; ep[1]=e[1]; ep[2]=e[2];
}

/* estimate F with the algebraic minimization algorithm from a set of matches in which outliers
 * are assumed to have been removed. F is assumed to contain an initial estimate (e.g. computed
 * with the eight point algorithm)
 *
 * Assumes that pts0, pts1, inliersidx & numInliers of indat have been set by the caller
 * but does not alter them in any way.
 */
static int estAlgMinEpipGeom(double F[NUM_FPARAMS], struct Fdata *indat, int normalize, int verbose)
{
register int i;
double (*inpts0)[2], (*inpts1)[2], L0[9], L1[9];
struct Fdatalg datalg;
register int *inliersidx;
int numInliers, ret;
double opts[LM_OPTS_SZ], info[LM_INFO_SZ], epip[3];

  if(indat->numInliers<NUM_FPARAMS) return FUNDEST_ERR;  // too few matches

  numInliers=indat->numInliers;
  datalg.cfs=(double **)malloc(numInliers*sizeof(double *));
  if(!datalg.cfs){
    fprintf(stderr, "Error: not enough memory for 'datalg.cfs' in estAlgMinEpipGeom()\n");
    exit(1);
  }

  /* one "big" malloc instead of several "small" ones */
  datalg.cfs[0]=(double *)malloc(numInliers*NUM_FPARAMS*sizeof(double)); /* one equation per point */
  if(!datalg.cfs[0]){
    fprintf(stderr, "Error: not enough memory for 'datalg.cfs[0]' in estAlgMinEpipGeom()\n");
    exit(1);
  }
  for(i=1; i<numInliers; ++i)
    datalg.cfs[i]=datalg.cfs[i-1] + NUM_FPARAMS;

  inliersidx=indat->inliersidx;
  inpts0=indat->pts0; inpts1=indat->pts1;

  if(normalize){
    double (*pts0)[2], (*pts1)[2];

    pts0=(double (*)[2])malloc(numInliers*sizeof(double[2]));
    pts1=(double (*)[2])malloc(numInliers*sizeof(double[2]));
    if(!pts0 || !pts1){
      fprintf(stderr, "Memory allocation request failed in estAlgMinEpipGeom()\n");
      exit(1);
    }

    /* copy inliers */
    for(i=0; i<numInliers; ++i){
      register int j;

      j=inliersidx[i];
      pts0[i][0]=inpts0[j][0]; pts0[i][1]=inpts0[j][1];
      pts1[i][0]=inpts1[j][0]; pts1[i][1]=inpts1[j][1];
    }
    fundest_normalize2DPts(pts0, pts0, numInliers, L0);
    fundest_normalize2DPts(pts1, pts1, numInliers, L1);
    fundest_normalizeF(F, L0, L1, F);

    for(i=0; i<numInliers; ++i)
      calcFundLinCoeffs2(pts0[i], pts1[i], datalg.cfs[i]);

    free(pts0);
    free(pts1);
  }
  else
    for(i=0; i<numInliers; ++i)
      calcFundLinCoeffs2(inpts0[inliersidx[i]], inpts1[inliersidx[i]], datalg.cfs[i]);

  datalg.prevepip[0]=datalg.prevepip[1]=datalg.prevepip[2]=0.0;

  /* minimize the algebraic error as function of the epipole */
  opts[0]=LM_INIT_MU; opts[1]=1E-12; opts[2]=1E-12; opts[3]=1E-15;
  opts[4]=LM_DIFF_DELTA;

  fundest_epipfromF(F, epip, NULL);
  // note zero measurement vector below
  ret=dlevmar_dif(algerror, epip, NULL, 3, numInliers, 1000, opts, info, NULL, NULL, (void *)&datalg); // no Jacobian

  fundest_min_AGx_norm(NULL, NULL, 0, 0, 0, NULL); // cleanup

  /* retrieve F */
  memcpy(F, datalg.fm, NUM_FPARAMS*sizeof(double));

  if(normalize)
    fundest_denormalizeF(F, L0, L1, F);

  inftynorm3x3(F); // ensure unit infty norm

  free(datalg.cfs[0]);
	free(datalg.cfs);

  return (ret>=0)? FUNDEST_OK : FUNDEST_ERR;
}
/**********************************************************************************************/

/* robust, non-linear fundamental matrix estimation from "nmatches" matched point features,
 * possibly including outliers. "inpts0", "inpts1" contain the matched point coordinates,
 * "inlPcent" is the expected percentage of inliers (>=0.5), "f" contains the estimated F
 * matrix upon return, "normalize" is a flag specifying whether input coordinates should
 * be normalized according to Hartley's suggestion, "howto" specifies whether a non-linear
 * refinement step shoud follow the linear estimation (8 point algorithm) and its cost
 * function (see fundest.h for appropriate values), "idxOutliers" points to sufficiently
 * large memory which upon return is set to the indices of the detected outlying points
 * (pass NULL if don't care), "nbOutliers" contains the number of outliers, "verbose"
 * specifies the verbosity level
 */
int fundest(double (*inpts0)[2], double (*inpts1)[2], int nmatches, double inlPcent, double F[NUM_FPARAMS],
            int normalize, int howto, int *idxOutliers, int *nbOutliers, int verbose)
{
register int i, j;
int isSqr=1, maxNbSol=1;
double gate=2.0, premResid=-1.0, sampleProb=0.99, outlierThresh;
int *outliersMap, ret, **sets=NULL, nbSets=0;
double (*pts0)[2], (*pts1)[2], L0[9], L1[9];
struct Fdata dat;

  if(howto==FUNDEST_NO_EST){
    fprintf(stderr, "\"FUNDEST_NO_EST\" option is meaningless in fundest()\n");
    return FUNDEST_ERR;
  }

  if(nmatches<NUM_LINFMATCHES) return FUNDEST_ERR;  // too few matches

  if(normalize){
    pts0=(double (*)[2])malloc(nmatches*sizeof(double[2]));
    pts1=(double (*)[2])malloc(nmatches*sizeof(double[2]));
    if(!pts0 || !pts1){
      fprintf(stderr, "Memory allocation request failed in fundest()\n");
      exit(1);
    }
    fundest_normalize2DPts(inpts0, pts0, nmatches, L0);
    fundest_normalize2DPts(inpts1, pts1, nmatches, L1);
  }
  else{
    pts0=inpts0;
    pts1=inpts1;
  }

  dat.cfs=(double **)malloc(nmatches*sizeof(double *));
  if(!dat.cfs){
    fprintf(stderr, "Error: not enough memory for 'dat.cfs' in fundest()\n");
    exit(1);
  }

  /* one "big" malloc instead of several "small" ones */
  dat.cfs[0]=(double *)malloc(nmatches*NUM_FPARAMS*sizeof(double)); /* one equation per point */
  if(!dat.cfs[0]){
    fprintf(stderr, "Error: not enough memory for 'dat.cfs[0]' in fundest()\n");
    exit(1);
  }
  for(i=1; i<nmatches; ++i)
    dat.cfs[i]=dat.cfs[i-1] + NUM_FPARAMS;

  for(i=0; i<nmatches; ++i)
    calcFundLinCoeffs2(pts0[i], pts1[i], dat.cfs[i]);

  nbSets=2*lqs_numtries(NUM_LINFMATCHES, inlPcent, sampleProb); // twice those really necessary
  sets=lqs_allocsets(NUM_LINFMATCHES, nbSets);

#ifdef USE_BUCKETS
  nbSets=fundest_genRandomSetsWithBuckets(pts0, NUM_LINFMATCHES, nmatches, nbSets, sets);
#else
  nbSets=fundest_genRandomSetsNoBuckets(NUM_LINFMATCHES, nmatches, nbSets, sets);
#endif /* USE_BUCKETS */

  dat.pts0=pts0; dat.pts1=pts1;
  dat.M=NULL; dat.Mrows=0;

  if(!(outliersMap=(int *)malloc(nmatches*sizeof(int)))){
    fprintf(stderr, "Error: not enough memory for 'outliersMap' in fundest()\n");
    exit(1);
  }
  verbose=verbose>1;

#if 1
  j=lqsfit(nmatches, NUM_LINFMATCHES, sets, nbSets, fmatrixLinResidualsAlg, estLinEpipGeom,
            isSqr, verbose, maxNbSol, gate, premResid, NUM_FPARAMS, inlPcent, (void *)&dat,
            F, NULL, outliersMap, nbOutliers, &outlierThresh);
#else
  //outlierThresh=ransac_getthresh(0.8, 1); // assume s=.8, symmetric distance involves 2 squared terms
  outlierThresh=sqrt(0.8*0.8*6.634896601); // assume s=.8, symmetric distance involves 2 squared terms
  j=ransacfit(nmatches, NUM_LINFMATCHES, sets, nbSets, fmatrixLinResidualsGeom, estLinEpipGeom,
            isSqr, verbose, maxNbSol, outlierThresh, 0, NUM_FPARAMS, inlPcent, (void *)&dat,
            F, NULL, outliersMap, nbOutliers);
#endif

  if(verbose){
    fprintf(stderr, "Outlier threshold: %g\n", outlierThresh);
    fprintf(stderr, "fundest(): LQS fit returned %d, %d outliers [out of %d]\n", j, *nbOutliers, nmatches);
  }

  if(sets) lqs_freesets(sets);

  if(j!=0){
    dat.inliersidx=(int *)malloc((nmatches - *nbOutliers)*sizeof(int));
    if(!dat.inliersidx){
      fprintf(stderr, "Error: not enough memory for 'dat.inliersidx' in fundest()\n");
      exit(1);
    }

    for(i=j=0; i<nmatches; ++i)
      if(!outliersMap[i]) dat.inliersidx[j++]=i;

    /* LS estimation on inliers */
    dat.numInliers=nmatches - *nbOutliers;
    estLinEpipGeom(F, dat.numInliers, dat.inliersidx, (void *)&dat);

    /* enforce the rank-2 constraint */
    mk2Rank3x3(F);

    /* expose outliers */
    if(idxOutliers!=NULL)
      for(i=j=0; i<nmatches; ++i)
        if(outliersMap[i]) idxOutliers[j++]=i;

#if 0
    if(verbose){
      fputs("Outliers: ", stderr);
      for(i=j=0; i<nmatches; ++i)
        if(outliersMap[i]) fprintf(stderr, "%d ", i);
      fputc('\n', stderr);
    }
#endif

    ret=FUNDEST_OK;

#if 0
    /* include the following code fragment to print the (unnormalized) matching 2D points found to be inlying */
    for(i=0; i<dat.numInliers; ++i){
      printf("%.4lf %.4lf  %.4lf %.4lf\n", inpts0[dat.inliersidx[i]][0], inpts0[dat.inliersidx[i]][1],
                                          inpts1[dat.inliersidx[i]][0], inpts1[dat.inliersidx[i]][1]);
    }
#endif

  }
  else{ /* LQS failed */
    memset(F, 0, NUM_FPARAMS*sizeof(double));
    *nbOutliers=nmatches;
    dat.numInliers=0; /* makes sure any refinements below are avoided */
    dat.inliersidx=NULL;
    ret=FUNDEST_ERR;
  }

  estLinEpipGeom(NULL, 0, NULL, (void *)&dat); /* release memory */

  /* the DLT estimate has now been computed. Time for any refinements: algebraic error or non-linear */

  if(howto==FUNDEST_ALGMIN && dat.numInliers>=NUM_FPARAMS)
    estAlgMinEpipGeom(F, &dat, !_DONORM, verbose); // algebraic error minimization

  if(normalize){
    free(pts0);
    free(pts1);
    fundest_denormalizeF(F, L0, L1, F);
  }

#if 0
    /* estimate before refinement */
    for(i=0; i<NUM_FPARAMS; ++i){
      if(!(i%3)) printf("\n");
      printf("%g ", F[i]);
    }
    fflush(stdout);
#endif

  if(dat.numInliers>=NUM_FNLPARAMS && howto!=FUNDEST_8PT && howto!=FUNDEST_ALGMIN){
    /* use the unnormalized points for the non-linear refinement */
    dat.pts0=inpts0; dat.pts1=inpts1;

    refineEpipGeom(F, howto, &dat, verbose);
  }

  inftynorm3x3(F); // ensure unit infty norm

  if(dat.inliersidx) free(dat.inliersidx);
  free(dat.cfs[0]);
	free(dat.cfs);
  free(outliersMap);

  return ret;
}

/* estimate F with the eight point algorithm (DLT) from a set of matches in which outliers
 * are assumed to have been removed.
 *
 * Assumes that pts0, pts1, inliersidx & numInliers of indat have been set by the caller
 * but does not alter them in any way.
 */
static int est8PtEpipGeom(double F[NUM_FPARAMS], struct Fdata *indat, int normalize, int verbose)
{
register int i;
double (*inpts0)[2], (*inpts1)[2], (*pts0)[2], (*pts1)[2], L0[9], L1[9];
struct Fdata dat;
register int *inliersidx;
int ret;

  dat.numInliers=indat->numInliers;
  if(dat.numInliers<NUM_LINFMATCHES) return FUNDEST_ERR;  // too few matches
  dat.inliersidx=NULL;

  inliersidx=indat->inliersidx;
  inpts0=indat->pts0; inpts1=indat->pts1;

  pts0=(double (*)[2])malloc(dat.numInliers*sizeof(double[2]));
  pts1=(double (*)[2])malloc(dat.numInliers*sizeof(double[2]));
  if(!pts0 || !pts1){
    fprintf(stderr, "Memory allocation request failed in est8PtEpipGeom()\n");
    exit(1);
  }

  /* copy inliers */
  for(i=0; i<dat.numInliers; ++i){
    register int j;

    j=inliersidx[i];
    pts0[i][0]=inpts0[j][0]; pts0[i][1]=inpts0[j][1];
    pts1[i][0]=inpts1[j][0]; pts1[i][1]=inpts1[j][1];
  }

  if(normalize){
    fundest_normalize2DPts(pts0, pts0, dat.numInliers, L0);
    fundest_normalize2DPts(pts1, pts1, dat.numInliers, L1);
  }

  dat.cfs=(double **)malloc(dat.numInliers*sizeof(double *));
  if(!dat.cfs){
    fprintf(stderr, "Error: not enough memory for 'dat.cfs' in est8PtEpipGeom()\n");
    exit(1);
  }

  /* one "big" malloc instead of several "small" ones */
  dat.cfs[0]=(double *)malloc(dat.numInliers*NUM_FPARAMS*sizeof(double)); /* one equation per point */
  if(!dat.cfs[0]){
    fprintf(stderr, "Error: not enough memory for 'dat.cfs[0]' in est8PtEpipGeom()\n");
    exit(1);
  }
  for(i=1; i<dat.numInliers; ++i)
    dat.cfs[i]=dat.cfs[i-1] + NUM_FPARAMS;

  for(i=0; i<dat.numInliers; ++i)
    calcFundLinCoeffs2(pts0[i], pts1[i], dat.cfs[i]);

  dat.pts0=pts0; dat.pts1=pts1;
  dat.M=NULL; dat.Mrows=0;

  /* LS estimation on all points */
  ret=estLinEpipGeom(F, dat.numInliers, NULL, (void *)&dat);

  estLinEpipGeom(NULL, 0, NULL, (void *)&dat); /* release memory */

  /* enforce the rank-2 constraint */
  mk2Rank3x3(F);

  if(normalize)
    fundest_denormalizeF(F, L0, L1, F);

  inftynorm3x3(F); // ensure unit infty norm

  free(pts0);
  free(pts1);

  free(dat.cfs[0]);
	free(dat.cfs);

  return (ret==1)? FUNDEST_OK : FUNDEST_ERR;
}

/* Similar to fundest but with an initial estimate provided in F.
 * Just linear estimation or non-linear refinement using the inliers determined with F; no LMedS/RANSAC.
 * If inlThresh>0, use it as absolute thresh, otherwise interpret |inlThresh| as a fraction to
 * determine thresh as a percentile. In the second case (inlThresh<0), if inlThresh<-1, its
 * absolute fractional part is used as a fraction and its integer part as an extra threshold
 * that is compared with that obtained with the fraction and the stricter of the two is retained.
 * In this manner, specifying e.g. -2.7 for inlThresh implies that a threshold equal to the
 * minimum of 2.0 and the 70th percentile of distances will be used.
 */
int fundest_wie(double (*pts0)[2], double (*pts1)[2], int nmatches, double inlThresh, double F[NUM_FPARAMS],
            int howto, int *idxOutliers, int *nbOutliers, int verbose)
{
register int i, j;
struct Fdata dat;

  memset(&dat, 0, sizeof(struct Fdata)); // clear everything

  dat.pts0=pts0; dat.pts1=pts1;
  dat.inliersidx=(int *)calloc(nmatches, sizeof(int)); // worst case, initialization necessary!
  if(!dat.inliersidx){
    fprintf(stderr, "Error: not enough memory for 'dat.inliersidx' in fundest_wie()\n");
    exit(1);
  }

  if(inlThresh>=0.0){
    double thresh, nom, denom, lin[3];
    register double *pt0, *pt1;

    thresh=inlThresh*inlThresh;
    for(i=j=0; i<nmatches; ++i){
      /* squared distance from the epipolar lines in both images */
      pt0=pts0[i]; pt1=pts1[i];

      /* image 1 */
      lin[0]=F[0]*pt0[0] + F[1]*pt0[1] + F[2];
      lin[1]=F[3]*pt0[0] + F[4]*pt0[1] + F[5];
      lin[2]=F[6]*pt0[0] + F[7]*pt0[1] + F[8];

      nom=(lin[0]*pt1[0] + lin[1]*pt1[1] + lin[2]);
      nom=SQR(nom);
      denom=(SQR(lin[0]) + SQR(lin[1]));
      if(nom>denom*thresh) continue; // do not bother checking image 0

      /* image 0 */
      lin[0]=F[0]*pt1[0]+ F[3]*pt1[1]+ F[6];
      lin[1]=F[1]*pt1[0]+ F[4]*pt1[1]+ F[7];
      lin[2]=F[2]*pt1[0]+ F[5]*pt1[1]+ F[8];

      nom=(lin[0]*pt0[0] + lin[1]*pt0[1] + lin[2]);
      nom=SQR(nom);
      denom=(SQR(lin[0]) + SQR(lin[1]));
      if(nom>denom*thresh) continue; 

      dat.inliersidx[j++]=i;
    }
    dat.numInliers=j;
  }
  else{
    double *distances, *distances_cpy;
    double thresh, thresh2=DBL_MAX, dist, lin[3], tmp;
    register double *pt0, *pt1;

    distances=(double *)malloc(2*nmatches*sizeof(double));
    if(!distances){
      fprintf(stderr, "Error: not enough memory for 'distances' in fundest_wie()\n");
      exit(1);
    }
    distances_cpy=distances+nmatches;

    for(i=0; i<nmatches; ++i){
      /* squared distance from the epipolar lines in both images */
      pt0=pts0[i]; pt1=pts1[i];

      /* image 1 */
      lin[0]=F[0]*pt0[0] + F[1]*pt0[1] + F[2];
      lin[1]=F[3]*pt0[0] + F[4]*pt0[1] + F[5];
      lin[2]=F[6]*pt0[0] + F[7]*pt0[1] + F[8];

      tmp=(lin[0]*pt1[0] + lin[1]*pt1[1] + lin[2]);
      dist=SQR(tmp)/(SQR(lin[0]) + SQR(lin[1]));

      /* image 0 */
      lin[0]=F[0]*pt1[0]+ F[3]*pt1[1]+ F[6];
      lin[1]=F[1]*pt1[0]+ F[4]*pt1[1]+ F[7];
      lin[2]=F[2]*pt1[0]+ F[5]*pt1[1]+ F[8];

      tmp=(lin[0]*pt0[0] + lin[1]*pt0[1] + lin[2]);
      dist+=SQR(tmp)/(SQR(lin[0]) + SQR(lin[1]));
      distances[i]=distances_cpy[i]=0.5*dist;
    }

    inlThresh=-inlThresh;
    if(inlThresh>1.0){
      /* inlThresh is not in [0, 1]. Extract its fractional part
       * to use as a fraction and combine the threshold computed
       * from it with its integer part
       */
      thresh2=floor(inlThresh);
      inlThresh=inlThresh-thresh2;
      thresh2*=thresh2;
    }

    thresh=lqs_kth_smallest(distances_cpy, nmatches, (int)(nmatches*inlThresh));

    if(thresh2<thresh) thresh=thresh2; // use the stricter of the two

    for(i=j=0; i<nmatches; ++i)
      if(distances[i]<=thresh) dat.inliersidx[j++]=i;

    dat.numInliers=j;
    free(distances);
  }

  verbose=verbose>1;
  switch(howto){
    case FUNDEST_NO_EST: // no estimation, inlier detection only
      break;

    case FUNDEST_8PT:
      est8PtEpipGeom(F, &dat, _DONORM, verbose); // F linearly estimated from the inliers detected with the provided F
      break;

    case FUNDEST_ALGMIN:
      estAlgMinEpipGeom(F, &dat, _DONORM, verbose); // algebraic error minimization for the inliers detected with the provided F
      break;

    case FUNDEST_EPIP_DIST:
    case FUNDEST_SAMPSON_ERROR:
      if(dat.numInliers>=NUM_FNLPARAMS)
        refineEpipGeom(F, howto, &dat, verbose); // non-linear refinement of provided F
      break;

    default:
      fprintf(stderr, "Unknown \'howto\' %d in fundest_wie()\n", howto);
      exit(1);
  }

  inftynorm3x3(F); // ensure unit infty norm

  /* expose outliers; assumes idxOutliers is sorted in ascending order */
  if(idxOutliers!=NULL){
    register int k;

    for(i=j=k=0; i<nmatches; ++i)
      if(dat.inliersidx[j]!=i) idxOutliers[k++]=i;
      else ++j;
  }

  free(dat.inliersidx);
  *nbOutliers=nmatches-dat.numInliers;

  return FUNDEST_OK;
}


/***** utility functions *****/

/* Compute the Root Mean Squared (RMS) and Root Median Squared (RMedS) reprojection error
 * pertaining to a camera matrix Pthat has been estimated from 2D-3D correspondences points.
 * Note that the RMS measure is sensitive to mismatched points while the RMedS is not
 */

#define _MEDIAN_(a, n) lqs_kth_smallest(a, n, (((n)&1)? ((n)/2) : (((n)/2)-1)))

void fundest_RMS_RMedS(double (*pts0)[2], double (*pts1)[2], int nmatches, double F[NUM_FPARAMS],
                 double *rms, double *rmeds)
{
register int i;
double *errors, sum, *pt0, *pt1, lin[3], dist, tmp;

  if((errors=(double *)malloc(nmatches*sizeof(double)))==NULL){
     fprintf(stderr, "Memory allocation request failed in fundest_RMS_RMedS()\n");
     exit(1);
  }

  for(i=0, sum=0.0; i<nmatches; ++i){
    /* symmetric distance from the epipolar lines in both images */
    pt0=pts0[i]; pt1=pts1[i];

    /* image 1 */
    lin[0]=F[0]*pt0[0] + F[1]*pt0[1] + F[2];
    lin[1]=F[3]*pt0[0] + F[4]*pt0[1] + F[5];
    lin[2]=F[6]*pt0[0] + F[7]*pt0[1] + F[8];

    tmp=(lin[0]*pt1[0] + lin[1]*pt1[1] + lin[2]);
    dist=SQR(tmp)/(SQR(lin[0]) + SQR(lin[1]));

    /* image 0 */
    lin[0]=F[0]*pt1[0]+ F[3]*pt1[1]+ F[6];
    lin[1]=F[1]*pt1[0]+ F[4]*pt1[1]+ F[7];
    lin[2]=F[2]*pt1[0]+ F[5]*pt1[1]+ F[8];

    tmp=(lin[0]*pt0[0] + lin[1]*pt0[1] + lin[2]);
    dist+=SQR(tmp)/(SQR(lin[0]) + SQR(lin[1]));
    dist*=0.5;

    sum+=errors[i]=dist;
  }

  *rms=sqrt(sum/(double)(nmatches));
  *rmeds=sqrt(_MEDIAN_(errors, nmatches));

  free(errors);
}

/* Compute the values dividing the reprojection errors into four equal parts.
 * The three values respectively cut off the lowest 25%, 50% and 75% of data.
 *
 */

void fundest_quartiles(double (*pts0)[2], double (*pts1)[2], int nmatches, double F[NUM_FPARAMS],
                      double *Q1, double *Q2, double *Q3)
{
register int i;
double *errors, *pt0, *pt1, lin[3], dist, tmp;

  if((errors=(double *)malloc(nmatches*sizeof(double)))==NULL){
     fprintf(stderr, "Memory allocation request failed in fundest_quartiles()\n");
     exit(1);
  }

  for(i=0; i<nmatches; ++i){
    /* symmetric distance from the epipolar lines in both images */
    pt0=pts0[i]; pt1=pts1[i];

    /* image 1 */
    lin[0]=F[0]*pt0[0] + F[1]*pt0[1] + F[2];
    lin[1]=F[3]*pt0[0] + F[4]*pt0[1] + F[5];
    lin[2]=F[6]*pt0[0] + F[7]*pt0[1] + F[8];

    tmp=(lin[0]*pt1[0] + lin[1]*pt1[1] + lin[2]);
    dist=SQR(tmp)/(SQR(lin[0]) + SQR(lin[1]));

    /* image 0 */
    lin[0]=F[0]*pt1[0]+ F[3]*pt1[1]+ F[6];
    lin[1]=F[1]*pt1[0]+ F[4]*pt1[1]+ F[7];
    lin[2]=F[2]*pt1[0]+ F[5]*pt1[1]+ F[8];

    tmp=(lin[0]*pt0[0] + lin[1]*pt0[1] + lin[2]);
    dist+=SQR(tmp)/(SQR(lin[0]) + SQR(lin[1]));
    dist*=0.5;

    errors[i]=dist;
  }

  if(Q1) *Q1=sqrt(lqs_kth_smallest(errors, nmatches, (nmatches>>2))); // 25%
  if(Q2) *Q2=sqrt(lqs_kth_smallest(errors, nmatches, (nmatches>>1))); // 50%
  if(Q3) *Q3=sqrt(lqs_kth_smallest(errors, nmatches, nmatches-(nmatches>>2))); // 75%

  free(errors);
}
