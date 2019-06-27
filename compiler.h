/////////////////////////////////////////////////////////////////////////////////
// 
//  Fundamental matrix estimation from 2D point matches
//  Copyright (C) 2002 - 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _COMPILER_H
#define _COMPILER_H

/* C compiler specific constants & macros */

/* lapack f77 routines name mangling */
#define F77_FUNC(func)    func ## _

/* note: intel's icc defines both __ICC & __INTEL_COMPILER.
 * Also, some compilers other than gcc define __GNUC__,
 * therefore gcc should be checked last
 */
#ifdef _MSC_VER
#define inline __inline // MSVC
#elif !defined(__ICC) && !defined(__INTEL_COMPILER) && !defined(__GNUC__)
#define inline // other than MSVC, ICC, GCC: define empty
#endif

#ifdef _MSC_VER
#define HYPOT2(a, b) _hypot((a), (b)) // MSVC
#elif defined(__ICC) || defined(__INTEL_COMPILER) || (defined(__GNUC__) && !defined(__MINGW32__)) // ICC, GCC: use hypot
#define HYPOT2(a, b) hypot((a), (b))
#else
#define HYPOT2(a, b) sqrt((a)*(a) + (b)*(b)) // other than MSVC, ICC, GCC: use sqrt(a^2 + b^2)
#endif

#ifdef _MSC_VER
#include <float.h>
#define POSEST_FINITE _finite // MSVC
#elif defined(__ICC) || defined(__INTEL_COMPILER) || defined(__GNUC__)
#define POSEST_FINITE finite // ICC, GCC
#else
#define POSEST_FINITE finite // other than MSVC, ICC, GCC, let's hope this will work
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#endif /* _COMPILER_H */
