/* ////////////////////////////////////////////////////////////////////////////////
// 
//  Matlab MEX file for fundest
//  Copyright (C) 2007-2012  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//////////////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <fundest.h>

#include "mex.h"

/**
#define DEBUG
**/

#define MAX(A, B)     ((A)>=(B)? (A) : (B))


/* display printf-style error messages in matlab */
static void matlabFmtdErrMsgTxt(char *fmt, ...)
{
char  buf[256];
va_list args;

	va_start(args, fmt);
	vsprintf(buf, fmt, args);
	va_end(args);

  mexErrMsgTxt(buf);
}

/* display printf-style warning messages in matlab */
static void matlabFmtdWarnMsgTxt(char *fmt, ...)
{
char  buf[256];
va_list args;

	va_start(args, fmt);
	vsprintf(buf, fmt, args);
	va_end(args);

  mexWarnMsgTxt(buf);
}

/* matlab matrices are in column-major, this routine converts them to row major for fundest */
static double *getTranspose(mxArray *Am)
{
int m, n;
double *At, *A;
register int i, j;

  m=mxGetM(Am);
  n=mxGetN(Am);
  A=mxGetPr(Am);
  At=mxMalloc(m*n*sizeof(double));

  for(i=0; i<m; i++)
    for(j=0; j<n; j++)
      At[i*n+j]=A[i+j*m];
  
  return At;
}

/*
[F, idxOutliers]=fundest(pts0, pts1, inlPcent, normalize, howto, verbose);
      pts0, pts1 are the input matching point pairs
      inlPcent is the expected fraction of inliers in the input points
      normalize specifies whether Hartley's normalization should be applied
                to input points prior to DLT estimation (default)
      howto controls how estimation is performed; can be one of:
               '8pt', 'algmin', 'epip_dist', 'sampson_err'
               or empty (implying '8pt', default)
               A synonym for '8pt' is 'norefine'
      verbose is optional

If an approximate estimate Fini of the fund. matrix is known, then fundest can be invoked as

[F, idxOutliers]=fundest(pts0, pts1, inlThresh, Fini, howto, verbose);
    inlThresh specifies how the outliers are to be determined; see fundest_wie() for details
    Remaining arguments are as above
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *Prhs[])
{
register int i;
int normalize, wie=0, howto, verbose, nbOutliers, *idxOutliers;
double (*pts0)[2], (*pts1)[2], inlPcent, F[NUM_FPARAMS];
register double *pdbl;
mxArray **prhs=(mxArray **)&Prhs[0];
int nmatches, len, status;

  /* parse input args; start by checking their number */
  if(nrhs<4 && nrhs>6)
    matlabFmtdErrMsgTxt("fundest: between 4 and 6 input arguments required (got %d).", nrhs);
  if(nlhs>2)
    matlabFmtdErrMsgTxt("fundest: at most 2 output arguments returned (got %d).", nlhs);
    
  /** pts0 **/
  /* first argument must be a two-column matrix */
  if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetN(prhs[0])!=2)
    matlabFmtdErrMsgTxt("fundest: first argument must be a two-column matrix (got %dx%d).", mxGetM(prhs[0]), mxGetN(prhs[0]));
  pts0=(double (*)[2])getTranspose(prhs[0]);
  nmatches=mxGetM(prhs[0]);

  /** pts1 **/
  /* second argument must be a two-column matrix */
  if(!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1])!=2)
    matlabFmtdErrMsgTxt("fundest: second argument must be a two-column matrix (got %dx%d).", mxGetM(prhs[1]), mxGetN(prhs[1]));
  if(mxGetM(prhs[1])!=nmatches)
    matlabFmtdErrMsgTxt("fundest: two first arguments should have the same number of rows (got %d and %d).", nmatches, mxGetM(prhs[1]));
  pts1=(double (*)[2])getTranspose(prhs[1]);

  /** inlPcent **/
  /* the third argument must be a scalar */
  if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
    mexErrMsgTxt("fundest: inlPcent must be a scalar.");
  inlPcent=mxGetScalar(prhs[2]);

  /** normalize or Fini **/
  /* check whether the fourth argument is a scalar or a 3x3 matrix */
  if(nrhs>=4 && (mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3]) && mxGetM(prhs[3])==1 && mxGetN(prhs[3])==1)){
    normalize=mxGetScalar(prhs[3])!=0.0;

    ++prhs;
    --nrhs;
  } else
  if(nrhs>=4 && (mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3]) && mxGetM(prhs[3])==3 && mxGetN(prhs[3])==3)){
    wie=1;
    pdbl=mxGetPr(prhs[3]);
    /* get transposed, i.e. row major */
    F[0]=pdbl[0]; F[1]=pdbl[3]; F[2]=pdbl[6];
    F[3]=pdbl[1]; F[4]=pdbl[4]; F[5]=pdbl[7];
    F[6]=pdbl[2]; F[7]=pdbl[5]; F[8]=pdbl[8];

    ++prhs;
    --nrhs;
  }
  else
    normalize=1;

  /** howto **/
  /* check whether fifth argument is a string */
  if(nrhs>=4 && mxIsChar(prhs[3])==1 && mxGetM(prhs[3])==1){
    char *str;

    /* examine supplied name */
    len=mxGetN(prhs[3])+1;
    str=mxCalloc(len, sizeof(char));
    status=mxGetString(prhs[3], str, len);
    if(status!=0)
      mexErrMsgTxt("fundest: not enough space. String is truncated.");

    for(i=0; str[i]; ++i)
      str[i]=tolower(str[i]);

    if(!strcmp(str, "8pt") || !strcmp(str, "norefine")) howto=FUNDEST_8PT;
    else if(!strcmp(str, "algmin")) howto=FUNDEST_ALGMIN;
    else if(!strcmp(str, "sampson_err")) howto=FUNDEST_SAMPSON_ERROR;
    else if(!strcmp(str, "epip_dist")) howto=FUNDEST_EPIP_DIST;
    else matlabFmtdErrMsgTxt("fundest: unknown minimization type '%s'.", str);

    mxFree(str);

    ++prhs;
    --nrhs;
  }
  else
    howto=FUNDEST_8PT;

  /** verbose **/
  /* the sixth argument must be a scalar */
  if(nrhs>=4){
    if(!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=1)
      mexErrMsgTxt("fundest: verbose must be a scalar.");
    verbose=mxGetScalar(prhs[3])!=0.0;
  }
  else
    verbose=0;

  if(nlhs>1) /* outlier indices should be returned */
    idxOutliers=mxMalloc(nmatches*sizeof(int));
  else
    idxOutliers=NULL;

  /* invoke fundest */
  if(!wie)
    fundest(pts0, pts1, nmatches, inlPcent, F, normalize, howto, idxOutliers, &nbOutliers, verbose);
  else // inlPcent is actually inlThresh...
    fundest_wie(pts0, pts1, nmatches, inlPcent, F, howto, idxOutliers, &nbOutliers, verbose);

  /* copy back returned results */
  /** F **/
  plhs[0]=mxCreateDoubleMatrix(3, 3, mxREAL);
  pdbl=mxGetPr(plhs[0]);
  /* return transposed, i.e. column major */
  pdbl[0]=F[0]; pdbl[1]=F[3]; pdbl[2]=F[6];
  pdbl[3]=F[1]; pdbl[4]=F[4]; pdbl[5]=F[7];
  pdbl[6]=F[2]; pdbl[7]=F[5]; pdbl[8]=F[8];

  /** idxOutliers **/
  if(nlhs>1){
    plhs[1]=mxCreateDoubleMatrix(nbOutliers, 1, mxREAL);
    pdbl=mxGetPr(plhs[1]);
    for(i=0; i<nbOutliers; ++i)
      *pdbl++=idxOutliers[i]+1; // convert from 0 to 1-based

    mxFree(idxOutliers);
  }

  /* cleanup */
  mxFree(pts0);
  mxFree(pts1);
}
