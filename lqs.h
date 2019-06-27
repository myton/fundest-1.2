/////////////////////////////////////////////////////////////////////////////////
// 
//  Fundamental matrix estimation from 2D point matches
//  Copyright (C) 2002 - 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _LQS_H
#define _LQS_H

extern int lqs_numtries(int p, double i, double P);
extern int **lqs_allocsets(int sizeSet, int nbSets);
extern void lqs_freesets(int **sets);
extern void lqs_initnoise();
extern void lqs_setrepeatablerandom(int flag);
extern double lqs_uniformnoise(double a, double b);

extern int lqsfit(int nbData, int sizeSet, int **sets, int nbSets,
                  void (*residual)(double *x, int nb, void *adata, double *resid),
                  int (*estimator)(double *x, int nb, int *indexes, void *adata),
                  int isResidualSqr, int verbose, int maxNbSol, double gate,
                  double prematureResidual, int dim, double percentageOfGoodData, void *adata,
                  double *estimate, int *bestSet, int *outliersMap, int *nbOutliers, double *outlierThresh);

extern double lqs_kth_smallest(double *a, int n, int k);

#endif /* _LQS_H */
