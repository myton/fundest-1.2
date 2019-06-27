/////////////////////////////////////////////////////////////////////////////////
// 
//  Geometric Robust Information Criterion (GRIC)
//  Copyright (C) 2011 - 2012  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FUNDEST_GRIC_FUND     0
#define FUNDEST_GRIC_AFUND    1
#define FUNDEST_GRIC_HOMO     2

/* Torr's geometric robust information criterion (GRIC)
 * for relative quality of fit assessment / model selection.
 *
 * "res" holds the squared residuals including outliers (i.e.,
 * the approximate **squared** geometric distances from
 * 4d joint-space points to the fitted manifold),
 * "sigma" is the assumed variance of the error, e.g. 0.8,
 * "n" the number of residuals and "model" specifies the
 * type of the fitted manifold.
 */
static double calc_GRIC(double *res, double sigma, int n, int model)
{
register int i;
int K, D; // number of parameters, dimension of the manifold
const int R=4; // data dimension (image point pairs)
const double sigmasq1=1.0/(sigma*sigma);
double lam3RD;
register double sum;

  switch(model){
  case 0: // fundamental 
    K=7; 
    D=3; 
    break;

  case 1: // affine fundamental
    K=4; 
    D=3;
      fprintf(stderr, "affine fundamental matrix not yet implemented in calc_GRIC()\n");
      exit(1);
    break;

  case 2: // homography
    K=8;
    D=2;
    break;

  default:
    fprintf(stderr, "unkown model %d specified to calc_GRIC()\n", model);
    exit(1);
  }

  lam3RD=2.0*(R-D); // lam3==2.0
  for(i=0, sum=0.0; i<n; ++i){
    double tmp;

    tmp=res[i]*sigmasq1;
    sum+=(tmp<=lam3RD)? tmp : lam3RD;
  }

  sum+=n*D*log((double)R) + K*log((double)(R*n));
	
  return sum;
}

/* squared geometric distance to the F manifold.
 *
 * see Fusiello's F_sampson_distance_sqr
 */
static void sampsonF_dsqr(double F[9], double (*pts0)[2], double (*pts1)[2], int npts, double *res)
{
register int i;
double *m0, *m1, Fm0[3], Ftm1[3], m1Fm0;

  for(i=0; i<npts; ++i){
    m0=pts0[i]; m1=pts1[i];

    /* F*m0 */
    Fm0[0]=F[0]*m0[0] + F[1]*m0[1] + F[2];
    Fm0[1]=F[3]*m0[0] + F[4]*m0[1] + F[5];
    Fm0[2]=F[6]*m0[0] + F[7]*m0[1] + F[8];

    /* F'*m1 */
    Ftm1[0]=F[0]*m1[0] + F[3]*m1[1] + F[6];
    Ftm1[1]=F[1]*m1[0] + F[4]*m1[1] + F[7];
    Ftm1[2]=F[2]*m1[0] + F[5]*m1[1] + F[8];

    /* m1'*F*m0 */
    m1Fm0=Fm0[0]*m1[0] + Fm0[1]*m1[1] +  Fm0[2];

    res[i]=m1Fm0*m1Fm0/(Fm0[0]*Fm0[0] + Fm0[1]*Fm0[1] + Ftm1[0]*Ftm1[0] + Ftm1[1]*Ftm1[1]);
  }
}

/* squared geometric distance to the H manifold.
 *
 * see vgg_H_sampson_distance_sqr.m, vgg_H_algebraic_distance.m
 */
static void sampsonH_dsqr(double H[9], double (*pts0)[2], double (*pts1)[2], int npts, double *res)
{
register int i;
double *m0, *m1;
double G0[3], G1[3], magG0, magG1, magG0G1, alpha, D1, D2, alg[2];

  for(i=0; i<npts; ++i){
    m0=pts0[i]; m1=pts1[i];

    G0[0]= H[0] - m1[0] * H[6];
    G0[1]= H[1] - m1[0] * H[7];
    G0[2]=-m0[0] * H[6] - m0[1] * H[7] - H[8];

    G1[0]= H[3] - m1[1] * H[6];
    G1[1]= H[4] - m1[1] * H[7];
    G1[2]=-m0[0] * H[6] - m0[1] * H[7] - H[8];

    magG0=sqrt(G0[0]*G0[0] + G0[1]*G0[1] + G0[2]*G0[2]);
    magG1=sqrt(G1[0]*G1[0] + G1[1]*G1[1] + G1[2]*G1[2]);
    magG0G1=G0[0]*G1[0] + G0[1]*G1[1];

    alpha=acos(magG0G1 /(magG0*magG1));

    /* algebraic distance */
    alg[0]=   m0[0]*H[0] + m0[1]*H[1] + H[2] -
       m1[0]*(m0[0]*H[6] + m0[1]*H[7] + H[8]);

    alg[1]=   m0[0]*H[3] + m0[1]*H[4] + H[5] -
       m1[1]*(m0[0]*H[6] + m0[1]*H[7] + H[8]);

    D1=alg[0]/magG0;
    D2=alg[1]/magG1;

    res[i]=(D1*D1 + D2*D2 - 2.0*D1*D2*cos(alpha))/sin(alpha);
  }
}

/* Given a set of matching points (including outliers) and estimates of fundamental
 * and homography matrices F & H resp., computes Torr's GRIC
 *
 * A model with lower GRIC is more likely
 */
void fundest_GRIC(double (*pts0)[2], double (*pts1)[2], int nmatches,
                  double F[9], double H[9], double sigma, double *gricF, double *gricH)
{
double *res;

  res=(double *)malloc(nmatches*sizeof(double));
  if(!res){
    fprintf(stderr, "memory allocation request failed in fundest_GRIC()\n");
    exit(1);
  }

  if(gricF){
    sampsonF_dsqr(F, pts0, pts1, nmatches, res);
    *gricF=calc_GRIC(res, sigma, nmatches, FUNDEST_GRIC_FUND);
  }

  if(gricH){
    sampsonH_dsqr(H, pts0, pts1, nmatches, res);
    *gricH=calc_GRIC(res, sigma, nmatches, FUNDEST_GRIC_HOMO);
  }

  free(res);
}
