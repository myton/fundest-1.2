/////////////////////////////////////////////////////////////////////////////////
// 
//  Fundamental matrix estimation from 2D point matches
//  Copyright (C) 2002 - 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _UTIL_H
#define _UTIL_H

/* linalg.c */
extern int fundest_min_Ax_normEIG(double *A, int m, int n, double *x);
extern int fundest_min_Ax_normSVD(double *A, int m, int n, double *x);
extern int fundest_min_AGx_norm(double *A, double *G, int m, int n, int r, double *x);
extern int fundest_matinvLU(double *A, double *B, int m);

/* norm.c */
extern void fundest_normalize2DPts(double (*pts)[2], double (*npts)[2], int numpts, double T[9]);
extern void fundest_normalizeF(double F[9], double T1[9], double T2[9], double nF[9]);
extern void fundest_denormalizeF(double nF[9], double T1[9], double T2[9], double F[9]);

/* buckets.c */
extern int fundest_genRandomSetsNoBuckets(int sizeSet, int nbData, int nbSets, int **sets);
extern int fundest_genRandomSetsWithBuckets(double (*pts)[2], int sizeSet, int nbData, int nbSets, int **sets);

#endif /* _UTIL_H */
