/*
/////////////////////////////////////////////////////////////////////////////////
// 
//  Fundamental matrix estimation from 2D point matches
//  Copyright (C) 2002 - 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////
*/

#ifndef _FUNDEST_H
#define _FUNDEST_H

/* define the following if you want to build a DLL with MSVC */
/**
#define DLL_BUILD 
**/

#ifdef __cplusplus
extern "C" {
#endif

#define FUNDEST_VERSION    "1.2 (Feb. 2016)"

/* estimation options & non-linear refinement cost functions */
#define FUNDEST_NO_EST           0 /* no estimation; fundest_wie() only */
#define FUNDEST_8PT              1 /* 8-point linear estimation */
#define FUNDEST_ALGMIN           2 /* algebraic minimization */
#define FUNDEST_EPIP_DIST        3 /* non-linear refinement using distance from epipolar lines */
#define FUNDEST_SAMPSON_ERROR    4 /* non-linear refinement using Sampson error */

#define FUNDEST_NO_NLN_REFINE    FUNDEST_8PT /* obsolete, for backwards compatibility only */

#define NUM_FPARAMS              9 /* #params involved in linear estimation */

/* use as: extern FUNDEST_API_MOD int FUNDEST_CALL_CONV func(...) */
#if defined(DLL_BUILD) && defined(_MSC_VER) /* build DLLs with MSVC only! */
#define FUNDEST_API_MOD    __declspec(dllexport)
#define FUNDEST_CALL_CONV  __cdecl
#else /* define empty */
#define FUNDEST_API_MOD 
#define FUNDEST_CALL_CONV
#endif /* DLL_BUILD && _MSC_VER */

#define FUNDEST_ERR     -1
#define FUNDEST_OK       0

/* fundest.c */
extern int fundest(double (*pts0)[2], double (*pts1)[2], int nmatches, double inlPcent,
                 double *F, int normalize, int howto, int *idxOutliers, int *nbOutliers, int verbose);

extern int fundest_wie(double (*pts0)[2], double (*pts1)[2], int nmatches, double inlThresh, double F[NUM_FPARAMS],
                int howto, int *idxOutliers, int *nbOutliers, int verbose);

/* gric.c */
extern void fundest_GRIC(double (*pts0)[2], double (*pts1)[2], int nmatches,
                      double F[9], double H[9], double sigma, double *gricF, double *gricH);

/* misc.c */
extern void fundest_epipfromF(double F[9], double *e, double *ep);
extern void fundest_EfromFK(double F[9], double K1[9], double K2[9], double E[9]);
extern void fundest_RtfromE(double E[9], double R1[9], double t1[3], double R2[9], double t2[3]);

/* RMS & RMedS errors for F */
extern void fundest_RMS_RMedS(double (*pts0)[2], double (*pts1)[2], int nmatches, double F[NUM_FPARAMS],
                             double *rms, double *rmeds);

/* first, second and third quartiles for F */
void fundest_quartiles(double (*pts0)[2], double (*pts1)[2], int nmatches, double F[NUM_FPARAMS],
                  double *Q1, double *Q2, double *Q3);

#ifdef __cplusplus
}
#endif

#endif /* _FUNDEST_H */
