/////////////////////////////////////////////////////////////////////////////////
// 
//  Fundamental matrix estimation from 2D point matches
//  Copyright (C) 2002 - 2015  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "compiler.h"

/* SVD */
extern int F77_FUNC(dgesvd)(char *jobu, char *jobvt, int *m, int *n,
                            double *a, int *lda, double *s, double *u, int *ldu,
                            double *vt, int *ldvt, double *work, int *lwork, int *info);

/* Extract the epipoles from the fundamental matrix.
 * e is the epipole in the 1st image and ep that in the 2nd,
 * i.e. F*e=0 and ep'*F=0
 * Any of them can be NULL, in which case no assignment is
 * made
 */
void fundest_epipfromF(double F[9], double *e, double *ep)
{
double Ft[9];
double S[3], U[9], Vt[9];

  /* transpose F */
  Ft[0]=F[0]; Ft[1]=F[3]; Ft[2]=F[6];
  Ft[3]=F[1]; Ft[4]=F[4]; Ft[5]=F[7];
  Ft[6]=F[2]; Ft[7]=F[5]; Ft[8]=F[8];

  /* find left/right null values (epipoles) */
#if 0 // use custom 3x3 SVD
  svd3(U, S, Vt, Ft);
  /* svd3 actually returns V; for compatibility with LAPACK which returns V^T,
   * the computed Vt is transposed in place to yield the true Vt
   */
  val=Vt[1]; Vt[1]=Vt[3]; Vt[3]=val;
  val=Vt[2]; Vt[2]=Vt[6]; Vt[6]=val;
  val=Vt[5]; Vt[5]=Vt[7]; Vt[7]=val;
#else // use generic LAPACK SVD
  {
  int n=3, info, lwork=32; /* probably too much here...*/
  double work[32];

  F77_FUNC(dgesvd)("S", "S", &n, &n, Ft, &n, S, U, &n, Vt, &n, work, &lwork, &info);

  if(info<0){
    fprintf(stderr, "LAPACK error: illegal value for argument %d of dgesvd in fundest_epipfromF()\n", -info);
    exit(1);
  }
  else if(info>0){
    fprintf(stderr, "LAPACK error: dbdsqr did not converge in fundest_epipfromF();\n%d %s", info,
      "superdiagonals of an intermediate bidiagonal form did not converge to zero\n");
    exit(1);
  }
  }
#endif

  /* note that dgesvd sorts singular values, thus guarantees
   * that the smallest (i.e. null) one is S[2]
   */

  /* left null vector */
  if(ep){
    ep[0]=U[6]; ep[1]=U[7]; ep[2]=U[8];
  }

  /* right null vector */
  if(e){
    e[0]=Vt[2]; e[1]=Vt[5]; e[2]=Vt[8];
  }
}

/* essential matrix from F matrix and calibration information.
 * E=K2'*F*K1
 *
 * If K2!=NULL, K2 is assumed equal to K1 (E=K1'*F*K1)
 */
void fundest_EfromFK(double F[9], double K1[9], double K2[9], double E[9])
{
register int i;
double tmp[9];

  if(K2==NULL) K2=K1;

  /* tmp=F*K1 */
  for(i=0; i<3; ++i){
    tmp[i*3  ]=F[i*3]*K1[0] + F[i*3+1]*K1[3  ] + F[i*3+2]*K1[2*3  ];
    tmp[i*3+1]=F[i*3]*K1[1] + F[i*3+1]*K1[3+1] + F[i*3+2]*K1[2*3+1];
    tmp[i*3+2]=F[i*3]*K1[2] + F[i*3+1]*K1[3+2] + F[i*3+2]*K1[2*3+2];
  }
  
  /* E=K2'*tmp */
  for(i=0; i<3; ++i){
    E[i*3  ]=K2[i]*tmp[0] + K2[1*3+i]*tmp[1*3  ] + K2[2*3+i]*tmp[2*3  ];
    E[i*3+1]=K2[i]*tmp[1] + K2[1*3+i]*tmp[1*3+1] + K2[2*3+i]*tmp[2*3+1];
    E[i*3+2]=K2[i]*tmp[2] + K2[1*3+i]*tmp[1*3+2] + K2[2*3+i]*tmp[2*3+2];
  }
}

/* decompose an essential matrix E to a combination of rotation and translation.
 *
 * E is assumed to be a proper essential matrix. i.e. has two equal singular
 * values and the third is zero. Since t can be determined up to scale, its
 * sign is unknown. Therefore, there are two possible solutions for t, each
 * associated with a rotation matrix. Thus, for the given E, there are two
 * R,t pairs, i.e. R1,t1 and R2,t2. Since -E is also a valid esential
 * matrix, R2,t1 and R1,t2 are also valid pairs, for a total of 4 possible
 * solutions. Of these, only one is physically plausible, i.e. gives rise
 * to reconstructed points which are all in front of both cameras.
 *
 * Computation is SVD-free, as described by Horn in "Recovering Baseline
 * and Orientation from 'Essential' Matrix".
 */
void fundest_RtfromE(double E[9], double R1[9], double t1[3], double R2[9], double t2[3])
{
double EEt[9], trace, s;
double Cf[9], TxE[9];
int i;

  /* EEt=E*E' */
  EEt[0]=E[0]*E[0] + E[1]*E[1] + E[2]*E[2];
  EEt[1]=E[0]*E[3] + E[1]*E[4] + E[2]*E[5];
  EEt[2]=E[0]*E[6] + E[1]*E[7] + E[2]*E[8];
  EEt[3]=E[3]*E[0] + E[4]*E[1] + E[5]*E[2];
  EEt[4]=E[3]*E[3] + E[4]*E[4] + E[5]*E[5];
  EEt[5]=E[3]*E[6] + E[4]*E[7] + E[5]*E[8];
  EEt[6]=E[6]*E[0] + E[7]*E[1] + E[8]*E[2];
  EEt[7]=E[6]*E[3] + E[7]*E[4] + E[8]*E[5];
  EEt[8]=E[6]*E[6] + E[7]*E[7] + E[8]*E[8];

  /* EEt=1/2*trace(E*E')*eye(3)-E*E' */
  trace=0.5*(EEt[0]+EEt[4]+EEt[8]);
  EEt[0]=trace-EEt[0];
  EEt[1]=-EEt[1];
  EEt[2]=-EEt[2];
  EEt[3]=-EEt[3];
  EEt[4]=trace-EEt[4];
  EEt[5]=-EEt[5];
  EEt[6]=-EEt[6];
  EEt[7]=-EEt[7];
  EEt[8]=trace-EEt[8];

  i=(EEt[0]>EEt[4])? ((EEt[0]>EEt[8])? 0 : 2) : ((EEt[4]>EEt[8])? 1 : 2); // max row
  s=sqrt(EEt[i*(3+1)]);
  t1[0]=EEt[i*3]/s; t1[1]=EEt[i*3+1]/s; t1[2]=EEt[i*3+2]/s;
  s=1.0/(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);

  /* Cf=cofactors(E) */
  Cf[0]= E[4]*E[8]-E[7]*E[5];
  Cf[1]=-E[8]*E[3]+E[6]*E[5];
  Cf[2]= E[7]*E[3]-E[6]*E[4];
  Cf[3]=-E[1]*E[8]+E[7]*E[2];
  Cf[4]= E[0]*E[8]-E[2]*E[6];
  Cf[5]=-E[0]*E[7]+E[1]*E[6];
  Cf[6]= E[1]*E[5]-E[4]*E[2];
  Cf[7]=-E[0]*E[5]+E[3]*E[2];
  Cf[8]= E[0]*E[4]-E[3]*E[1];

  /* TxE=skewsym(t1)*E */
  TxE[0]=-t1[2]*E[3]+t1[1]*E[6];
  TxE[1]=-t1[2]*E[4]+t1[1]*E[7];
  TxE[2]=-t1[2]*E[5]+t1[1]*E[8];
  TxE[3]= t1[2]*E[0]-t1[0]*E[6];
  TxE[4]= t1[2]*E[1]-t1[0]*E[7];
  TxE[5]= t1[2]*E[2]-t1[0]*E[8];
  TxE[6]=-t1[1]*E[0]+t1[0]*E[3];
  TxE[7]=-t1[1]*E[1]+t1[0]*E[4];
  TxE[8]=-t1[1]*E[2]+t1[0]*E[5];

  /* R1=(Cofactors(E) - skewsym(t1)*E)./(t1'*t1) */
  R1[0]=(Cf[0]-TxE[0])*s;
  R1[1]=(Cf[1]-TxE[1])*s;
  R1[2]=(Cf[2]-TxE[2])*s;
  R1[3]=(Cf[3]-TxE[3])*s;
  R1[4]=(Cf[4]-TxE[4])*s;
  R1[5]=(Cf[5]-TxE[5])*s;
  R1[6]=(Cf[6]-TxE[6])*s;
  R1[7]=(Cf[7]-TxE[7])*s;
  R1[8]=(Cf[8]-TxE[8])*s;

  /* R2=(Cofactors(E) + skewsym(t1)*E)./(t1'*t1) */
  R2[0]=(Cf[0]+TxE[0])*s;
  R2[1]=(Cf[1]+TxE[1])*s;
  R2[2]=(Cf[2]+TxE[2])*s;
  R2[3]=(Cf[3]+TxE[3])*s;
  R2[4]=(Cf[4]+TxE[4])*s;
  R2[5]=(Cf[5]+TxE[5])*s;
  R2[6]=(Cf[6]+TxE[6])*s;
  R2[7]=(Cf[7]+TxE[7])*s;
  R2[8]=(Cf[8]+TxE[8])*s;

  /* normalize translation */
  s=sqrt(s);
  t1[0]*=s; t1[1]*=s; t1[2]*=s;

  t2[0]=-t1[0]; t2[1]=-t1[1]; t2[2]=-t1[2];
}

