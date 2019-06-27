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
#include <math.h>
#include <float.h>

#include "util.h"
#include "lqs.h"

#undef USE_UNIQUE_SETS

//#define USE_UNIQUE_SETS

#define NUM_X_BUCKETS 2
#define NUM_Y_BUCKETS 2

#define MAX_TRIES       60

struct bucket{
  int nummatches;	            /* number of matches in bucket */
  int *indxs;                 /* index into points array */
};

struct selection{
  struct bucket *buck;        /* points to corresponding bucket */
  int selflag;	              /* indicates whether the corresponding bucket has been selected */
  double interval_width;	    /* width of the corresponding interval in [0, 1] */
};

static struct bucket buckets[NUM_X_BUCKETS][NUM_Y_BUCKETS];

static void shellsort_int(int *v, int n)
{
register int gap, i, j, jg;
register int tmp;

  for(gap=n>>1; gap>0; gap>>=1) // gap=n/2, gap/=2
    for(i=gap; i<n; ++i)
      for(j=i-gap; j>=0 && v[j]>v[jg=(j+gap)]; j-=gap){
        tmp=v[j];
        v[j]=v[jg];
        v[jg]=tmp;
      }
}

/* Sets generation without buckets */

/*
 * generate one subset making sure that all its elements are distinct
 *
 */
static void genOneSubsetNoBuckets(double h, int sizeSet, int subset[])
{
register int i, j, r;

  for(i=0; i<sizeSet; ++i){
resample:
    r=subset[i]=(int)lqs_uniformnoise(0.0, h);
    for(j=0; j<i; ++j)
      if(r==subset[j]) goto resample; /* element already in subset, try another one */
  }
}

/*
 * generate random subsets without using buckets
 *
 * returns the number of random sets actually generated
 */
int fundest_genRandomSetsNoBuckets(int sizeSet, int nbData, int nbSets, int **sets)
{
  register int i, j, k;
  int *work, *subset;

#ifdef USE_UNIQUE_SETS
  int in, out, tries=(3*nbSets)>>1; /* 3*nbSets/2, maximum number of tries */
#endif /* USE_UNIQUE_SETS */

  lqs_initnoise();

  if((work=(int *)malloc(nbData*sizeof(int)))==NULL){
    fprintf(stderr, "Memory allocation request failed in fundest_genRandomSetsNoBuckets()!\n");
    exit(1);
  }
  /* init work with "inside-out" Fisher-Yates shuffle */
  for(i=0; i<nbData; ++i){
    j=(int)lqs_uniformnoise(0, i); // random int in {0 ... i}
    work[i]=work[j];
    work[j]=i;
  }

  for(i=0; i<nbSets;  ){
      /* generate a subset */
      subset=sets[i];
#if 1
      for(j=0; j<sizeSet; ++j){
        k=(int)lqs_uniformnoise(j, nbData-0.01);  // select index from {j ... nbData-1} ...
        subset[j]=work[k]; work[k]=work[j]; work[j]=subset[j]; // ...and swap it with work[j]
      }
#else
      /* following code adapted from genOneSubsetNoBuckets() which is identical to lqs_gensubset() */
      for(j=0; j<sizeSet; ++j){
        register int r;

resample:
          r=subset[j]=(int)lqs_uniformnoise(0.0, nbData-0.5);
          for(k=0; k<j; ++k)
	          if(r==subset[k]) goto resample; /* element already in subset, try another one */
      }
#endif

#ifdef USE_UNIQUE_SETS
      shellsort_int(subset, sizeSet);
      /* verify whether the subset is already retained in sets */
      for(j=in=0; j<i && in==0; ++j){
	      for(k=out=0; k<sizeSet && out==0; ++k)
	        if(subset[k]!=sets[j][k]) out=1;
	      if(out==0) in=1;
	    }
      if(in==0) ++i;

      if(--tries==0){ /* premature exit, too many tries failed */
        //fprintf (stderr,"Warning: premature return in fundest_genRandomSetsNoBuckets(), too many tries failed!\n");
        break;
	    }
#else
      ++i;
#endif /* USE_UNIQUE_SETS */
  }

  free(work);

  return i;
}

/* Sets generation with buckets */

/* For more details on buckets, see INRIA RR-2273, pp. 19 */

/* initializes buckets and returns 1 if enough nonempty buckets exist, 0 otherwise */
static int initBuckets(double (*pts)[2], int npts, int sizeSet)
{
register int i, bx, by;
register struct bucket *bptr;
int numnonempty;
double coord, maxu, minu, maxv, minv, bucket_width, bucket_height;

	/* compute min-max corner coordinates */
	maxu=maxv=DBL_MIN;
	minu=minv=DBL_MAX;

	for(i=0; i<npts; ++i){
    coord=pts[i][0];
    if(coord > maxu) maxu=coord;
    if(coord < minu) minu=coord;

    coord=pts[i][1];
    if(coord > maxv) maxv=coord;
    if(coord < minv) minv=coord;
	}

  maxu++; maxv++;
  minu--; minv--;

	/* init bucket array */
	for(bx=0; bx<NUM_X_BUCKETS; ++bx)
		for(by=0; by<NUM_Y_BUCKETS; ++by){
			buckets[bx][by].nummatches=0;
			if(!(buckets[bx][by].indxs=(int *)malloc(npts*sizeof(int)))){
        fprintf(stderr, "Memory allocation request failed in initBuckets()!\n");
        exit(1);
      }
		}

	/* count points in buckets */
	bucket_width=(maxu-minu)/(double)(NUM_X_BUCKETS);
	bucket_height=(maxv-minv)/(double)(NUM_Y_BUCKETS);

	for(i=0; i<npts; ++i){
    coord=pts[i][0];
		bx=(int)((coord-minu)/bucket_width);

    coord=pts[i][1];
		by=(int)((coord-minv)/bucket_height);

    bptr=buckets[bx] + by;

    bptr->indxs[bptr->nummatches]=i;
    ++(bptr->nummatches);
    //printf("POINT %d goes to bucket %d %d\n", i, bx, by);
	}

	/* count nonempty buckets */
  for(i=numnonempty=0, bptr=buckets[0]; i<NUM_X_BUCKETS*NUM_Y_BUCKETS; ++i, ++bptr){
      if(bptr->nummatches==0) continue;

			++numnonempty;
  }

	if(numnonempty<sizeSet){
		fprintf(stderr, "Warning: Number of nonempty buckets in initBuckets() (%d) is smaller than required (%d)!\n%s\n",
						numnonempty, sizeSet, "Buckets will not be used; to enable them reduce their number");
		return 0;
	}

  return 1;
}

/* generate a random subset using buckets. No check is performed for 
 * avoiding cases where a point is contained more than once in the
 * computed subset
 */
static void genOneSetWithBuckets(double nbData, int sizeSet, int subset[])
{
register int i, j, tries;
register struct bucket *bptr;
int     nselected;
double 	anum, sum;
struct  selection selected_bucks[NUM_X_BUCKETS*NUM_Y_BUCKETS];

	/* prepare to start selecting buckets */
  for(i=j=0, bptr=buckets[0]; i<NUM_X_BUCKETS*NUM_Y_BUCKETS; ++i, ++bptr){
			if(bptr->nummatches==0) continue;

			selected_bucks[j].buck=bptr;
			selected_bucks[j].interval_width=((double)(bptr->nummatches))/nbData;
			selected_bucks[j].selflag=0;
			++j;
		}

	nselected=0;
	while(nselected<sizeSet){
		tries=MAX_TRIES;
    do{
		  anum=lqs_uniformnoise(0.0, 0.999);
      --tries;

		  i=-1; sum=0.0;
		  do{
			  ++i;
			  sum+=selected_bucks[i].interval_width;
		  } while(sum<anum);

    } while(selected_bucks[i].selflag && tries!=0);

    /* bucket `i' is either not selected or the `tries' limit has been reached */

		selected_bucks[i].selflag=1;

		/* choose a point from the selected bucket */
    j=(int)lqs_uniformnoise(0.0, (double)((selected_bucks[i].buck)->nummatches-0.5));
    subset[nselected]=(selected_bucks[i].buck)->indxs[j];
		++nselected;
	}
}

/*
 * generate random subsets using buckets
 *
 * returns the number of random sets actually generated
 */
int fundest_genRandomSetsWithBuckets(double (*pts)[2], int sizeSet, int nbData, int nbSets, int **sets)
{
  register int i, j;
  int *subset, ret;
  void (*func)(double, int, int*);

#ifdef USE_UNIQUE_SETS
  register int k;
  int in, out, tries=(3*nbSets)>>1; /* 3*nbSets/2, maximum number of tries */
#endif /* USE_UNIQUE_SETS */

  lqs_initnoise();

  if(initBuckets(pts, nbData, sizeSet))
    func=genOneSetWithBuckets;
  else
    func=genOneSubsetNoBuckets;
  
  for(i=0; i<nbSets;  ){
      /* generate a subset */
      subset=sets[i];
      (*func)(nbData, sizeSet, subset);

#ifdef USE_UNIQUE_SETS
      shellsort_int(subset, sizeSet);
      /* verify whether the subset is already contained in sets */
      for(j=in=0; j<i && in==0; ++j){
	      for(k=out=0; k<sizeSet && out==0; ++k)
	        if(subset[k]!=sets[j][k]) out=1;
	      if(out==0) in=1;
	    }
      if(in==0) ++i;

      if(--tries==0){ /* premature exit, too many tries failed */
        //fprintf (stderr,"Warning: premature return in fundest_genRandomSetsWithBuckets(), too many tries failed!\n");
	      break;
	    }
#else
      ++i;
#endif /* USE_UNIQUE_SETS */
  }
  ret=i;

  /* free bucket indices */
	for(i=0; i<NUM_X_BUCKETS; ++i)
		for(j=0; j<NUM_Y_BUCKETS; ++j)
      free(buckets[i][j].indxs);

  return ret;
}
