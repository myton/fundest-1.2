/////////////////////////////////////////////////////////////////////////////////
// 
//  Fundamental matrix estimation from 2D point matches
//  Copyright (C) 2002 - 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _CALC_FUND_COEFFS_H
#define _CALC_FUND_COEFFS_H

extern void calcFundLinCoeffs2(double m1[2], double m2[2], double eq[9]);

extern void calcEpipDist_12(double m1[2], double m2[2], double f8[8], double d[1]);
extern void calcEpipDistJac_12(double m1[2], double m2[2], double f8[8], double d_grad[8]);

extern void calcEpipDist_13(double m1[2], double m2[2], double f8[8], double d[1]);
extern void calcEpipDistJac_13(double m1[2], double m2[2], double f8[8], double d_grad[8]);

extern void calcEpipDist_23(double m1[2], double m2[2], double f8[8], double d[1]);
extern void calcEpipDistJac_23(double m1[2], double m2[2], double f8[8], double d_grad[8]);

extern void calcEpipConstrSampsonDist_12(double m1[2],double m2[2],double f8[8], double err[1]);
extern void calcEpipConstrSampsonDistGrads_12(double m1[2],double m2[2],double f8[8], double grads[8]);

extern void calcEpipConstrSampsonDist_13(double m1[2],double m2[2],double f8[8], double err[1]);
extern void calcEpipConstrSampsonDistGrads_13(double m1[2],double m2[2],double f8[8], double grads[8]);

extern void calcEpipConstrSampsonDist_23(double m1[2],double m2[2],double f8[8], double err[1]);
extern void calcEpipConstrSampsonDistGrads_23(double m1[2],double m2[2],double f8[8], double grads[8]);

#endif /* _CALC_FUND_COEFFS_H */
